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

/**
 * When a type is repeat (and not the first compressed split)
 * make sure to compute the stats over all the repeated blocks.
 */
typedef struct {
    U32 end;      /* once sequence past the end */
    char modes[4]; /* lit, ml, ll, of: basic, compressed, rle, repeat, any */
} blockSplit_t;

typedef struct {
    U32 idx;
    symbolEncodingType_e mode;
} pred_t;

typedef struct {
  U16 freqs[256]; /* maxNbSeq < 2^16-1 */
  double accum;   /* accumulator */
  double maxCost; /* Maximum cost allowed for this window */
  U32 endIdx;        /* one index past the end */
  U32 endSeq;        /* one sequence past the end */
} window_t;


/**
 * Returns #splits or an error code
 */
size_t ZSTD_blockSplit(ZSTD_CCtx const* zc, blockSplit_t* splits);

#endif /* BLOCK_SPLITTER_H */
