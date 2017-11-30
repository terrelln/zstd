#ifndef BLOCK_SPLITTER_H
#define BLOCK_SPLITTER_H

#include "zstd_internal.h"

/**
 * maxNbSplits: Maximum number of splits to check. If it is < #seq, then we will
 * skip every N sequences for for the smallest N that makes the #seq that we
 * check less than maxNbSplits.
 *
 * nbWindows: # of windows to
 */


typedef enum {
    st_lit = 0, st_off, st_ml, st_ll, st_end
} splitType_e;

/**
 * When a type is repeat (and not the first compressed split)
 * make sure to compute the stats over all the repeated blocks.
 */
typedef struct {
    U32 end;      /* Number of sequences in the block */
    char modes[4]; /* lit, ml, ll, of: basic, compressed, rle, repeat, any */
} blockSplit_t;

typedef struct {
    U32 idx;
    symbolEncodingType_e mode;
} pred_t;

typedef struct {
  U16 freqs[256];    /* maxNbSeq < 2^16-1 */
  U32 accum;         /* accumulator */
  U32 maxCost;       /* Maximum cost allowed for this window */
  U32 endIdx;        /* one index past the end */
  U32 endSeq;        /* one sequence past the end */
} window_t;

typedef struct {
  window_t* windows;
  blockSplit_t* splits;
  pred_t* pred;
  U32* offsets;
  U32* minCost;
  size_t nbWindows;
  size_t maxNbSplits;
} blockSplitState_t;

/* Returns the amount of space we need to allocate */
size_t ZSTD_sizeof_blockSplitState(size_t maxNbSeq);

/* Takes a pointer, initializes the state, and returns tail pointer */
void* ZSTD_blockSplitState_init(blockSplitState_t *state, void *ptr,
                                size_t maxNbSeq);

/**
 * Returns #splits or an error code
 */
size_t ZSTD_blockSplit(ZSTD_CCtx* zc);

size_t ZSTD_getNbSeqInChunk(blockSplit_t const* splits, size_t const nbSplits,
                            splitType_e type);
size_t ZSTD_getNbSeqInBlock(blockSplit_t const* split);
symbolEncodingType_e ZSTD_getBlockMode(blockSplit_t const* split,
                                       splitType_e type);

#endif /* BLOCK_SPLITTER_H */
