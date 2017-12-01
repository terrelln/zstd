#include "block_splitter.h"
#include "zstd_compress.h"
#include "zstd_internal.h"
#include <stdlib.h>
#include <math.h>

static pred_t const kEmptyPred = {-1, -1};

static U32 ZSTD_log2_U64(U64 const input)
{
  if (input == 0)
    return 0;
  return 63 - __builtin_clzll(input);
}

static U32 ZSTD_log2_U32(U32 const input)
{
  if (input == 0)
    return 0;
  return 31 - __builtin_clz(input);
}

#define ZSTD_log2(input) (sizeof(input) == 4 ? ZSTD_log2_U32(input) : assert(sizeof(input) == 8), ZSTD_log2_U64(input))


/**
 * Translate a sequence index into an offset into the symbols array.
 */
static size_t computeOffset(size_t const seq, U32 const *const offsets,
                            int const type)
{
    switch (type) {
    case st_lit:
        /* Literals can only be split on sequence boundaries.
         * Each sequence can have 0 or more literals.
         * Index is the sequence index, so it needs to be translated into an
         * offset into the literals buffer, which is symbols.
         */
        assert(offsets != NULL);
        return offsets[seq];
    case st_off:
    case st_ml:
    case st_ll:
        /* There is exactly one offset, match length, and literal length per
         * sequence, so index can used directly.
         */
        assert(offsets == NULL);
        return seq;
    default:
        return assert(0), 0;
    }
}

/**
 * Extends the window to index.
 * Precondition: window->end < index && index < # symbols.
 *
 * Templated over literals, offsets, match length, and literal lengths
 */
static void extendWindow(window_t* const window, size_t const idx,
                         size_t const seq, BYTE const* const symbols,
                         U32 const *const offsets, int const type)
{
    BYTE const *symbol = symbols + computeOffset(window->endSeq, offsets, type);
    BYTE const *const symbolEnd = symbols + computeOffset(seq, offsets, type);

    assert(window->endIdx < idx);
    assert(window->endSeq < seq);
    assert(symbol <= symbolEnd);
    /* We maintain window->accum = the sum of F[s] * log2(F[s]) for each symbol.
     * To update it we have to subtract the old frequency of the updated symbol,
     * and add back in the new frequency.
     */
    for (; symbol < symbolEnd; ++symbol) {
        U16 *const freq = &window->freqs[*symbol];

        window->accum -= *freq * ZSTD_log2(*freq);
        ++*freq;
        window->accum += *freq * ZSTD_log2(*freq);
    }
    window->endIdx = idx;
    window->endSeq = seq;
}

static void popWindows(window_t* const windows, size_t const nbWindows,
                       size_t const seq, BYTE const* const symbols,
                       U32 const *const offsets, int const type)
{
    BYTE const *const symbolsBegin =
        symbols + computeOffset(seq, offsets, type);
    BYTE const *const symbolsEnd =
        symbols + computeOffset(seq + 1, offsets, type);
    size_t i;

    assert(symbolsBegin <= symbolsEnd);

    for (i = 0; i < nbWindows; ++i) {
        BYTE const* symbol;

        assert(seq < windows[i].endSeq);
        if (seq + 1 >= windows[i].endSeq) {
            /* The window is empty, so we can clear the state */
            memset(windows[i].freqs, 0, sizeof(windows[i].freqs));
            windows[i].accum = 0;
            windows[i].endSeq = seq + 1;
        }
        for (symbol = symbolsBegin; symbol < symbolsEnd; ++symbol) {
            U16 *const freq = &windows[i].freqs[*symbol];

            windows[i].accum -= *freq * ZSTD_log2(*freq);
            --*freq;
            windows[i].accum += *freq * ZSTD_log2(*freq);
        }
    }
}

/**
 * Get the symbols array from the seqStore_t.
 */
static BYTE const* getSymbols(seqStore_t const* const seqStore, int const type)
{
    switch (type) {
    case st_lit: return seqStore->litStart;
    case st_off: return seqStore->ofCode;
    case st_ml:  return seqStore->mlCode;
    case st_ll:  return seqStore->llCode;
    default:     return assert(0), NULL;
    }
}

typedef struct {
    U32 cost;
    symbolEncodingType_e type;  /* basic, compressed, rle, repeat, any */
} cost_t;

/**
 * Compute the cost of a window.
 */
static cost_t computeCost(window_t const* window, size_t const seq,
                          U32 const *const offsets, size_t blockSwitchCost,
                          int const type) {
    size_t const windowBegin = computeOffset(seq, offsets, type);
    size_t const windowEnd = computeOffset(window->endSeq, offsets, type);
    size_t const windowSize = windowEnd - windowBegin;
    cost_t cost;
    assert(seq < window->endSeq);
    assert(windowEnd >= windowBegin);

    cost.cost = blockSwitchCost +
                windowSize * ZSTD_log2(windowSize) - window->accum;
    cost.type = set_compressed;

    assert(cost.cost > 0);
    DEBUGLOG(7, "windowSize = %u\twindowCost = %u\n", (U32)windowSize,
             cost.cost);
    return cost;
}

static size_t computeStepSize(size_t const nbSeq, size_t const maxNbSplits)
{
    return nbSeq / maxNbSplits + (nbSeq != maxNbSplits);
}

static double computeEpsilon(size_t const nbWindows,
                             size_t const blockSwitchCost, size_t const maxCost)
{
    return pow(maxCost / blockSwitchCost, 1.0 / nbWindows) - 1;
}

static blockSplit_t constructBlockSplit(size_t const endSeq,
                                        symbolEncodingType_e const mode,
                                        int const type)
{
    blockSplit_t result = {0, {set_repeat, set_repeat, set_repeat, set_repeat}};

    result.end = endSeq;
    result.modes[type] = mode;
    return result;
}

static size_t computeBlockCostAdjustment(blockSplit_t const* const splits, size_t const nbSplits, size_t const seq)
{
    size_t const kMinBlockCost = 6;
    size_t split;
    for (split = 0; split < nbSplits; ++split) {
        if (splits[split].end == seq)
            return kMinBlockCost;
    }
    return 0;
}

/* Block splitting is a single source shortest path (SSSP) problem on a DAG.
 * Data of length N has a node for each position, labeled 0..N. Each node
 * represents a possible split point. Each node i has directed edges to every
 * node greater than it, i+1, ..., N. A particular block split is a path from
 * node 0 to node N. Each edge (i, j) (j > i) represents a block from i to j.
 * Naturally, the cost of edge (i, j) is the cost to encode the block it
 * represents. We see that the solution to the SSSP problem from nodes 0 to N
 * is the optimal block splitting. We can solve the SSSP problem on a DAG in
 * O(V + E), assuming edge cost computations are constant. However, there are 2
 * problems:
 *
 * 1. There are O(N^2) edges, so the algorithm is at least O(N^2).
 * 2. The cost of (i, j), C(i, j), takes O(j - i) time to compute.
 *
 * Thus, we need to prune some edges, and approximate the edge cost, instead of
 * computing it exactly.
 *
 * Edge Pruning:
 *
 * We can prune some edges while maintaining a (1 + eps)-approximation of the
 * optimal solution. There is a minimum fixed cost for creating a new block
 * (e.g. the block header, the tables), which we will call F. We also know that
 * no block can cost more than 128 KB + 3 to encode, which we will call M.
 * For each cost in the series F, F(1 + eps), F(1 + eps)^2, ..., M, we will
 * maintain the first edge that costs no less than the series cost. Call the
 * maximum number of edges per node W. This guarantees a (1 + eps)-approximation
 * of the optimal solution while pruning the graph down to O(N/eps*log(M-F))
 * edges.
 *
 * Cost Approximation:
 *
 * We need to be able to answer the cost of edge (i, j), C(i, j), in O(1) time.
 * We are maintining W edges per node in the pruned graph. We will maintain W
 * windows, which we use to maintain state for each edge out of the current
 * node. Each window will be associated with a cost, and will expand until its
 * edge reaches that cost. When we move to the next node, we pop the earliest
 * node out of the window, and expand until its edge reaches its cost again.
 * We design the window such that it can answer C(i, j) in constant time, and
 * pushing the next node, and popping the earliest node are both O(1). Thus, we
 * can maintain the W windows, and add only a factor of W to the running time.
 * The final running time will be O(N * (log(M-F)/eps)^2).
 */
/* Template emulated for literals */
/* literals have offsets non-null */
FORCE_INLINE_TEMPLATE
size_t ZSTD_blockSplit_internal(seqStore_t const* const seqStore,
                                blockSplitState_t* const bs,
                                blockSplit_t* const splits,
                                size_t const maxNbSplits,
                                size_t const nbExistingSplits,
                                size_t const blockSwitchCost,
                                size_t const maxCost,
                                U32 const* const offsets,
                                int const type) {
    size_t const nbSeq = seqStore->sequences - seqStore->sequencesStart;
    size_t const stepSize = computeStepSize(nbSeq, maxNbSplits);
    size_t const nbSplits = nbSeq / stepSize;
    size_t const nbWindows = bs->nbWindows;
    window_t *const windows = bs->windows;
    pred_t *const pred = bs->pred;
    U32 *const minCost = bs->minCost;
    BYTE const *const symbols = getSymbols(seqStore, type);
    size_t idx;
    size_t seq;

    memset(windows, 0, nbWindows * sizeof(*windows));
    {
        double const epsilon =
            computeEpsilon(nbWindows, blockSwitchCost, maxCost);
        double windowMaxCost = blockSwitchCost;
        for (idx = 0; idx < nbWindows; ++idx) {
            windows[idx].maxCost = (U32)windowMaxCost;
            windowMaxCost *= 1 + epsilon;
        }
        assert(windowMaxCost >= maxCost - 1);
    }

    // Initialize minCost
    minCost[0] = 0.0;
    for (idx = 1; idx < nbSplits; ++idx) {
        minCost[idx] = (U32)-1;
    }
    if (kDebugLevel) {
      memset(pred, -1, nbSplits *sizeof(blockSplit_t));
    }
    /**
     * idx is the index into the DP arrays
     * seq is the sequences
     */
    // TODO: Figure out end codition for seq not divisible by stepSize.
    for (idx = 0, seq = 0; idx < nbSplits; ++idx, seq += stepSize) {
        size_t const blockCostAdjustment =
            computeBlockCostAdjustment(bs->splits, nbExistingSplits, seq);
        size_t lastIdxEnd = idx + 1;
        size_t lastSeqEnd = seq + 1;
        size_t windowIndex;

        assert(minCost[idx] != (U32)-1);
        assert(seq < nbSeq);

        for (windowIndex = 0; windowIndex < nbWindows; ++windowIndex) {
            window_t *const window = windows + windowIndex;
            int stop = 0;

            assert(lastSeqEnd < nbSeq);

            /* This window can be extended to the end of the previous window by
             * monotonicity.
             */
            if (window->endIdx < lastIdxEnd) {
                assert(window->endSeq < lastSeqEnd);
                extendWindow(window, lastIdxEnd, lastSeqEnd, symbols, offsets,
                             type);
                assert(window->endSeq == window->endIdx * stepSize);
            }

            for (;;) {
                cost_t const cost =
                    computeCost(window, seq, offsets,
                                blockSwitchCost - blockCostAdjustment, type);
                /* If this split is the minimum cost to get to window->endIdx
                 * so far save it.
                 */
                if (minCost[idx] + cost.cost < minCost[window->endIdx]) {
                    pred_t p = {idx, cost.type};

                    minCost[window->endIdx] = minCost[idx] + cost.cost;
                    pred[window->endIdx] = p;
                    DEBUGLOG(7, "pred[%u] = %u\n", (U32)window->endIdx,
                             (U32)idx);
                }
                /* If the window extends to the end of the block stop.
                 * The rest of the windows will have the same end, so we don't
                 * need to check them.
                 */
                if (window->endIdx == nbSplits) {
                    assert(window->endSeq > nbSeq - stepSize);
                    DEBUGLOG(7, "stop window %u\n", (U32)windowIndex);
                    stop = 1;
                    break;
                }
                /* If the window reached its maximum cost stop. */
                if (cost.cost > window->maxCost)
                    break;
                DEBUGLOG(7, "windows[%u].cost = %u\n", (U32)windowIndex,
                         cost.cost);
                /* Extend this window by one idx (stepSize seqs). */
                extendWindow(window, window->endIdx + 1,
                             window->endSeq + stepSize, symbols, offsets, type);
                assert(window->endSeq == window->endIdx * stepSize);
            }

            DEBUGLOG(
                6, "windowIndex %u\twindowSeqSize %u\twindowCost "
                   "%u\twindowMaxCost %u\n",
                (U32)windowIndex, (U32)(window->endSeq - seq),
                computeCost(window, seq, offsets, blockSwitchCost, type).cost,
                window->maxCost);
            /* If a window reached its maximum cost then skip the rest. */
            if (stop)
                break;
            lastIdxEnd = window->endIdx;
            lastSeqEnd = window->endSeq;
        }
        /* Remove the current sequence from every window. */
        popWindows(windows, nbWindows, seq, symbols, offsets, type);
    }

    {
        size_t nbPartitions = 1;
        size_t partition;
        /* Compute the length of the best path */
        for (idx = nbSplits - 1; idx != 0; idx = pred[idx].idx) {
            assert(idx < nbSplits);
            assert(minCost[idx] != (U32)-1);
            assert(pred[idx].idx != (U32)-1 &&
                   pred[idx].mode != (symbolEncodingType_e)-1);
            ++nbPartitions;
        }
        /* Fill the path */
        partition = nbPartitions;
        for (idx = nbSplits - 1; idx != 0; idx = pred[idx].idx) {
            assert(partition > 0);
            splits[--partition] =
                constructBlockSplit(idx * stepSize, pred[idx].mode, type);
        }
        assert(partition == 0);
        return nbPartitions;
    }
}

U32 const* initOffsets(U32* const offsets, seqDef const* const sequences,
                       size_t const nbSeq, size_t const litLength,
                       int const type)
{
    size_t seq;

    if (type != st_lit)
        return NULL;

    offsets[0] = 0;
    for (seq = 0; seq < nbSeq; ++seq) {
        offsets[seq + 1] = offsets[seq] + sequences[seq].litLength;
        assert(offsets[seq] <= offsets[seq + 1]);
    }
    offsets[seq + 1] = litLength;
    assert(offsets[seq] <= offsets[seq + 1]);

    return offsets;
}

/**
 * Computes the block switch cost in bits.
 */
static size_t computeBlockSwitchCost(int const type)
{
    size_t const kMaxBlockCost = 12;
    switch (type) {
    case st_lit: return kMaxBlockCost + ((size_t)80 << 3);
    case st_off: return kMaxBlockCost + (assert(0 && "TODO"), (size_t)80 << 3);
    case st_ml:  return kMaxBlockCost + (assert(0 && "TODO"), (size_t)80 << 3);
    case st_ll:  return kMaxBlockCost + (assert(0 && "TODO"), (size_t)80 << 3);
    default:     return assert(0), 0;
    }
}

/**
 * Computes the maximum cost of encoding the symbols of the given type.
 */
static size_t computeMaxCost(ZSTD_CCtx const* const zc, int const type)
{
    size_t const blockSwitchCost = computeBlockSwitchCost(type);
    size_t const nbSeq = zc->seqStore.sequences - zc->seqStore.sequencesStart;
    size_t const nbLits = zc->seqStore.lit - zc->seqStore.litStart;

    switch(type) {
    case st_lit: return blockSwitchCost + (nbLits << 3);
    // TODO: These bounds could probably be tighter because the range is
    //       smaller...
    case st_off: return blockSwitchCost + (nbSeq << 3);
    case st_ml:  return blockSwitchCost + (nbSeq << 3);
    case st_ll:  return blockSwitchCost + (nbSeq << 3);
    default:     return assert(0), 0;
    }
}

static int compareSplits(void const* lp, void const* rp)
{
    blockSplit_t const* const lhs = (blockSplit_t const*)lp;
    blockSplit_t const* const rhs = (blockSplit_t const*)rp;

    if (lhs->end < rhs->end)
        return -1;
    if (lhs->end > rhs->end)
        return 1;
    assert(lhs->end == rhs->end);
    return 0;
}

static int getType(blockSplit_t const split)
{
    int type;
    for (type = 0; type < 4; ++type)
        if (split.modes[type] != set_repeat)
            return type;
    return assert(0), 0;
}

static size_t mergeSplits(blockSplit_t* const splits, size_t const nbSplits)
{
    blockSplit_t* const endSplit = splits + nbSplits;
    blockSplit_t* outSplit = splits;
    blockSplit_t* inSplit = splits;
    /* Merge equal splits */
    while (++inSplit != endSplit) {
      if (outSplit->end == inSplit->end) {
          int const type = getType(*inSplit);
          assert(outSplit->modes[type] == set_repeat);
          outSplit->modes[type] = inSplit->modes[type];
      } else if (++outSplit != inSplit) {
        *outSplit = *inSplit;
      }
    }
    /* Merge repeats */
    return (outSplit - splits) + 1;
}

static blockSplit_t const* getFirstSplit(blockSplit_t const* splits,
                                         size_t const nbSplits,
                                         splitType_e type)
{
    blockSplit_t const* const end = splits + nbSplits;
    while (splits < end && splits->modes[type] == set_repeat)
        ++splits;
    return splits < end ? splits : NULL;
}

size_t ZSTD_getNbSeqInChunk(blockSplit_t const* splits, size_t const nbSplits,
                            splitType_e type)
{
    blockSplit_t const* const end = splits + nbSplits;
    size_t nbSeq = splits->end;
    while (++splits < end && splits->modes[type] == set_repeat) {
        nbSeq += splits->end;
    }
    return nbSeq;
}

size_t ZSTD_getNbSeqInBlock(blockSplit_t const* split)
{
    return split->end;
}

symbolEncodingType_e ZSTD_getBlockMode(blockSplit_t const* split,
                                       splitType_e type)
{

    return split->modes[type];
}

size_t ZSTD_blockSplit(ZSTD_CCtx* const zc)
{
    blockSplitState_t* const bs = &zc->blockSplitState;
    blockSplit_t* const splits = bs->splits;
    seqStore_t const* const seqStore = &zc->seqStore;
    seqDef const* const seqs = seqStore->sequencesStart;
    size_t const nbSeq = seqStore->sequences - seqStore->sequencesStart;
    size_t const litLength = seqStore->lit - seqStore->litStart;
    size_t nbSplits = 0;
    int type;

    // TODO: Maybe the order matters here, because we incorperate the existing
    // blocks. If literals get the biggest gain, maybe they should go last, or
    // first.
    for (type = st_lit; type <= st_lit; ++type) {
    // for (type = st_lit; type < st_end; ++type) {
        U32 const* const offsets = initOffsets(bs->offsets, seqs, nbSeq,
                                               litLength, type);
        size_t const blockSwitchCost = computeBlockSwitchCost(type);
        size_t const maxCost = computeMaxCost(zc, type);
        // TODO: Check that we don't overrun splits...
        nbSplits += ZSTD_blockSplit_internal(
            seqStore, bs, splits + nbSplits, bs->maxNbSplits, nbSplits,
            blockSwitchCost, maxCost, offsets, type);
    }
    // qsort(splits, nbSplits, sizeof(blockSplit_t), compareSplits);
    // nbSplits = mergeSplits(splits, nbSplits);
    for (type = st_lit + 1; type < st_end; ++type) {
        splits[0].modes[type] = set_any;
    }
    if (kDebugLevel) {
        size_t split;

        assert(splits[0].end > 0);
        for (split = 1; split < nbSplits; ++split) {
          assert(splits[split - 1].end < splits[split].end);
        }
    }
    // TODO: Attempt to join blocks that are very close together.
    //       To do that we need to find out the smallest legit block size we
    //       want to allow, and join blocks smaller than that.
    // Maybe not? If we can incorperate the cost into the algorithm...
    // Make splitting at block boundaries very attractive?
    return nbSplits;
}

size_t ZSTD_sizeof_blockSplitState(size_t maxNbSeq)
{
  blockSplitState_t state;
  return (size_t)ZSTD_blockSplitState_init(&state, NULL, maxNbSeq);
}

void* ZSTD_blockSplitState_init(blockSplitState_t* state, void* const ptr, size_t maxNbSeq)
{
  size_t const kNbWindows = 100;
  size_t const kMaxNbSplits = 10000;

  state->windows = (window_t*)ptr;
  state->splits = (blockSplit_t*)(state->windows + kNbWindows);
  state->pred = (pred_t*)(state->splits + kMaxNbSplits);
  state->offsets = (U32*)(state->pred + maxNbSeq);
  state->minCost = (U32*)(state->offsets + maxNbSeq);
  return state->minCost + maxNbSeq;
}

#if 0
double const kEpsilon = 0.3;
// TODO: Figure out average / p90 block switch costs.
// Check that is actually around 80 bytes, which experimentally works best.
double const kBlockSwitchCost = (size_t)80 << 3;
// TODO: Is it okay that this is 2^18, 2^20 is safe but maybe too large.
// This is too loose since most literals aren't actually 128 KB.
double const kMaxSinglePartitionCost = (size_t)1 << 20;

static double fastLog2(double const input) {
  return input == 0.0 ? 0.0 : log2(input);
}

static size_t computeNbWindows(size_t const litLength) {
  size_t nbWindows = 0;
  double const maxSinglePartitionCost = kBlockSwitchCost + (litLength << 3);
  double maxCost = kBlockSwitchCost;
  for (;;) {
    ++nbWindows;
    if (maxCost >= maxSinglePartitionCost) {
      break;
    }
    maxCost *= 1 + kEpsilon;
  }
  return nbWindows;
}

static void removeSplit(BYTE const *const literals, size_t const *const offsets,
                        window_t *const windows, size_t const nbWindows,
                        size_t const index) {
  BYTE const *const symbolsBegin = literals + offsets[index];
  BYTE const *const symbolsEnd = literals + offsets[index + 1];
  size_t i;

  for (i = 0; i < nbWindows; ++i) {
    BYTE const *symbol;
    // The window is empty
    if (index + 1 >= windows[i].end) {
      memset(windows[i].freqs, 0, sizeof(windows[i].freqs));
      windows[i].accum = 0;
      windows[i].end = index + 1;
      continue;
    }
    for (symbol = symbolsBegin; symbol < symbolsEnd; ++symbol) {
      U16 *const freq = windows[i].freqs + *symbol;

      windows[i].accum -= *freq * fastLog2(*freq);
      --*freq;
      windows[i].accum += *freq * fastLog2(*freq);
    }
  }
}

static void extendWindow(BYTE const *const literals,
                         size_t const *const offsets, size_t const nbSplits,
                         window_t *const window, size_t const index) {
  BYTE const *symbol = literals + offsets[window->end];
  BYTE const *const symbolsEnd = literals + offsets[index];

  assert(window->end < index);
  assert(index < nbSplits);
  assert(symbol < symbolsEnd);
  for (; symbol < symbolsEnd; ++symbol) {
    U16 *const freq = window->freqs + *symbol;

    window->accum -= *freq * fastLog2(*freq);
    ++*freq;
    window->accum += *freq * fastLog2(*freq);
  }
  window->end = index;
}

static size_t computeNbSplits(seqDef const *const sequences, size_t const nbSeq,
                              size_t const litLength) {
  size_t splits = 1;
  size_t i;
  size_t length = 0;
  for (i = 0; i < nbSeq; ++i) {
    if (sequences[i].litLength == 0) {
      continue;
    }
    length += sequences[i].litLength;
    ++splits;
  }
  if (length < litLength) {
    ++splits;
  }
  return splits;
}

double computeEpsilon(size_t const nbWindows, double const blockSwitchCost,
                      double const maxCost) {
    return pow(maxCost / blockSwitchCost, 1.0 / nbWindows) - 1;
}

FORCE_INLINE
size_t computeStepSize(size_t const nbSeq, size_t const maxNbSplits) {
    if (nbSeq > maxNbSplits) {
        size_t const stepSize = nbSeq / maxNbSplits + 1;
        assert((stepSize - 1) * maxNbSplits < nbSeq);
        return stepSize;
    }
    return 1;
}

/* Template emulated for literals */
/* literals have offsets non-null */
FORCE_INLINE
size_t ZSTD_blockSplit_internal(ZSTD_CCtx const *zc, blockSplit_t *splits,
                                BYTE const *const codes,
                                size_t const size,
                                double const blockSwitchCost,
                                double const maxCost,
                                size_t const *const offsets,
                                int type) {
    seqStore_t const *const seqStore = zc->seqStorePtr;
    size_t const nbSeq = seqStore->sequences - seqStore->sequencesStart;
    size_t const stepSize = computeNbSplits(nbSeq, zc->maxNbSplits);
    size_t const nbSplits = nbSeq / stepSize;
    size_t const nbWindows = zc->nbBlockSplitWindows;
    double const epsilon = computeEpsilon(nbWindows, blockSwitchCost, maxCost);
    window_t *const windows = zc->blockSplitWindows;
    blockSplit_t *const pred = zc->blockSplitPred;
    double *const minCost = zc->blockSplitMinCost;
    int i;
    int seq;

    memset(windows, 0, nbWindows * sizeof(*windows));
    {
        double windowMaxCost = blockSwitchCost;
        for (i = 0; i < nbWindows; ++i) {
            windows[i].maxCost = windowMaxCost;
            windowMaxCost *= 1 + epsilon;
        }
        assert(windowMaxCost >= maxCost - 1);
    }

    // Initialize minCost
    minCost[0] = 0.0;
    for (i = 1; i < nbSplits; ++i) {
        minCost[i] = INFINITY;
    }
#ifdef ZSTD_DEBUG
    // Initialize pred for debugging
    for (i = 0; i < nbSplits; ++i) {
        pred[i] = (size_t)-1;
    }
#endif
    for (i = 0, seq = 0; i < nbSplits; ++i, seq += stepSize) {
        assert(minCost[i] != INFINITY);
        assert(seq < nbSeq);
    }
}

size_t *ZSTD_block_split(seqStore_t const *seqStore, size_t *nbPartitions) {
  BYTE const *const literals = seqStore->litStart;
  size_t const litLength = seqStore->lit - literals;
  seqDef const *const sequences = seqStore->sequencesStart;
  size_t const nbSeq = seqStore->sequences - seqStore->sequencesStart;
  size_t const nbSplits = computeNbSplits(sequences, nbSeq, litLength);
  size_t const nbWindows = computeNbWindows(litLength);

  // TODO: Can I combine pred and offsets?
  size_t *const offsets = malloc(nbSplits * sizeof(size_t));
  size_t *const pred = malloc(nbSplits * sizeof(size_t));
  double *const minCost = malloc(nbSplits * sizeof(double));
  window_t *const windows = calloc(nbWindows, sizeof(window_t));
  size_t *partitions = NULL;
  size_t i;

  if (!offsets || !windows || !pred || !minCost) {
    goto out;
  }

  DEBUGLOG(2, "nbSeq: %zu\n", nbSeq);
  DEBUGLOG(2, "nbSplits: %zu\n", nbSplits);
  DEBUGLOG(2, "nbWindows: %zu\n", nbWindows);
  DEBUGLOG(2, "window fixed cost = %f\n", kBlockSwitchCost);

  // Set the fixed cost for a split for each window
  {
    double maxCost = kBlockSwitchCost;
    for (i = 0; i < nbWindows; ++i) {
      windows[i].accum = 0;
      windows[i].maxCost = maxCost;
      maxCost *= 1 + kEpsilon;
    }
    assert(maxCost >= kBlockSwitchCost + (litLength << 3));
  }

  // Compute all of the split points, including the end of the sequence
  {
    size_t offsetIdx = 1;
    offsets[0] = 0;
    for (i = 0; i < nbSeq; ++i) {
      if (sequences[i].litLength == 0) {
        continue;
      }
      offsets[offsetIdx] = offsets[offsetIdx - 1] + sequences[i].litLength;
      assert(offsets[offsetIdx - 1] < offsets[offsetIdx]);
      ++offsetIdx;
    }
    if (offsetIdx < nbSplits) {
      offsets[offsetIdx] = litLength;
      assert(offsets[offsetIdx - 1] < offsets[offsetIdx]);
      ++offsetIdx;
    }
    assert(offsetIdx == nbSplits);
    assert(offsets[nbSplits - 1] == litLength);
  }

  // Initialize minCost
  minCost[0] = 0.0;
  for (i = 1; i < nbSplits; ++i) {
    minCost[i] = INFINITY;
  }
#ifndef NDEBUG
  // Initialize pred for debugging
  for (i = 0; i < nbSplits; ++i) {
    pred[i] = (size_t)-1;
  }
#endif
  // Populate minCost from left to right
  for (i = 0; i < nbSplits - 1; ++i) {
    size_t lastEnd = i + 1;
    size_t windowIndex;
    assert(minCost[i] != INFINITY);
    DEBUGLOG(3, "i: %zu\n", i);

    for (windowIndex = 0; windowIndex < nbWindows; ++windowIndex) {
      window_t *const window = windows + windowIndex;
      int stop = 0;

      assert(lastEnd < nbSplits);

      if (window->end < lastEnd) {
        extendWindow(literals, offsets, nbSplits, window, lastEnd);
      }

      for (;;) {
        double const cost = computeCost(offsets, i, window);

        if (minCost[i] + cost < minCost[window->end]) {
          // Relax
          minCost[window->end] = minCost[i] + cost;
          pred[window->end] = i;
          DEBUGLOG(4, "pred[%zu] = %zu\n", window->end, i);
        }

        if (window->end == nbSplits - 1) {
          DEBUGLOG(4, "stop window %zu\n", windowIndex);
          stop = 1;
          break;
        }

        if (cost > window->maxCost) {
          // Can't grow this window anymore.
          break;
        }
        DEBUGLOG(4, "windows[%zu].cost = %f\n", windowIndex, cost);

        extendWindow(literals, offsets, nbSplits, window, window->end + 1);
      }

      DEBUGLOG(
          3,
          "windowIndex %zu\twindowSize %zu\twindowCost %f\twindowMaxCost %f\n",
          windowIndex, window->end - i, computeCost(offsets, i, window),
          window->maxCost);

      if (stop) {
        // No more windows can be extended this loop
        break;
      }
      lastEnd = window->end;
    }

    removeSplit(literals, offsets, windows, nbWindows, i);
  }

  // Compute the path length
  *nbPartitions = 1;
  {
    size_t pos = nbSplits - 1;
    while (pos != 0) {
      assert(pos < nbSplits);
      assert(minCost[pos] != INFINITY);
      assert(pred[pos] != (size_t)-1);
      ++*nbPartitions;
      pos = pred[pos];
    }
  }
  // Fill the path
  {
    size_t pos = nbSplits - 1;
    size_t *partition;
    partitions = malloc(*nbPartitions * sizeof(size_t));
    if (!partitions) {
      goto out;
    }
    partition = partitions + *nbPartitions;
    while (partition != partitions) {
      *--partition = offsets[pos];
      pos = pred[pos];
    }
  }
  --*nbPartitions;
  {
    // Log the path
    DEBUGLOG(2, "partitions of %zu:\n", seqStore->lit - literals);
    for (i = 1; i <= *nbPartitions; ++i) {
      DEBUGLOG(2, "\t%zu\n", partitions[i]);
    }
  }
// Translate the path to sequence numbers.
// {
//   size_t seq = 0;
//   size_t seqOffset = 0;
//   for (i = 0; i < *nbPartitions; ++i) {
//     DEBUGLOG(2, "partitions[i] = %zu\n", partitions[i]);
//     while (seqOffset < partitions[i]) {
//       if (seq == nbSeq) {
//         seqOffset = litLength;
//         assert(i == *nbPartitions - 1);
//         break;
//       }
//       seqOffset += sequences[seq].litLength;
//       ++seq;
//     }
//     assert(seqOffset == partitions[i]);
//     partitions[i] = seq;
//   }
//   DEBUGLOG(2, "seqOffset %zu\tlitLength%zu\n", seqOffset, litLength);
//   assert(seqOffset == litLength);
//   // DEBUGLOG(2, "seq %zu\tnbSeq%zu\n", seq, nbSeq);
//   // assert(seq == nbSeq);
// }
// {
//   // Log the path
//   DEBUGLOG(2, "partitions of %zu:\n", seqStore->lit - literals);
//   for (i = 0; i < *nbPartitions; ++i) {
//     DEBUGLOG(2, "\t%zu\n", partitions[i]);
//   }
// }
out:
  free(offsets);
  free(pred);
  free(minCost);
  free(windows);
  return partitions;
}
#endif
