/**
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

/*-**************************************
*  Tuning parameters
****************************************/
typedef double COVER_score_t;

/*-*************************************
*  Dependencies
***************************************/
#include <stdio.h>  /* fprintf */
#include <stdlib.h> /* malloc, free, qsort */
#include <string.h> /* memset */
#include <time.h>   /* clock */

#include "mem.h"           /* read */
#include "zstd_internal.h" /* includes zstd.h */
#ifndef ZDICT_STATIC_LINKING_ONLY
#define ZDICT_STATIC_LINKING_ONLY
#endif
#include "zdict.h"

/*-*************************************
*  Constants
***************************************/
#define COVER_MAX_SAMPLES_SIZE ((U32)-1)
#define COVER_ENTROPY_TABLES_SIZE 200

/*-*************************************
*  Console display
***************************************/
static int g_displayLevel = 2;
#define DISPLAY(...)                                                           \
  {                                                                            \
    fprintf(stderr, __VA_ARGS__);                                              \
    fflush(stderr);                                                            \
  }
#define DISPLAYLEVEL(l, ...)                                                   \
  if (g_displayLevel >= l) {                                                   \
    DISPLAY(__VA_ARGS__);                                                      \
  } /* 0 : no display;   1: errors;   2: default;  3: details;  4: debug */

#define DISPLAYUPDATE(l, ...)                                                  \
  if (g_displayLevel >= l) {                                                   \
    if ((clock() - g_time > refreshRate) || (g_displayLevel >= 4)) {           \
      g_time = clock();                                                        \
      DISPLAY(__VA_ARGS__);                                                    \
      if (g_displayLevel >= 4)                                                 \
        fflush(stdout);                                                        \
    }                                                                          \
  }
static const clock_t refreshRate = CLOCKS_PER_SEC * 15 / 100;
static clock_t g_time = 0;

/**
 * Returns the sum of the sample sizes.
 */
static size_t COVER_sum(const size_t *samplesSizes, unsigned nbSamples) {
  size_t sum = 0;
  unsigned i;
  for (i = 0; i < nbSamples; ++i) {
    sum += samplesSizes[i];
  }
  return sum;
}

typedef struct {
  U32 begin;
  U32 end;
  COVER_score_t score;
} COVER_segment_t;

typedef struct {
  const BYTE *samples;
  const size_t *offsets;
  size_t nbSamples;
  U32 *suffix;
  U32 *freqs;
  U32 *dmerAt;
  U32 *activeDmers;
  size_t size;
  COVER_params_t parameters;
} COVER_ctx_t;

/* We need a global context for qsort... */
static COVER_ctx_t *g_ctx = NULL;

/**
 * Returns -1 if the dmer at lp is less than the dmer at rp.
 * Return 0 if the dmers at lp and rp are equal.
 * Returns 1 if the dmer at lp is greater than the dmer at rp.
 */
static int COVER_cmp(COVER_ctx_t *ctx, const void *lp, const void *rp) {
  const U32 lhs = *(const U32 *)lp;
  const U32 rhs = *(const U32 *)rp;
  return memcmp(ctx->samples + lhs, ctx->samples + rhs, ctx->parameters.d);
}

/**
 * Same as COVER_cmp() except ties are broken by pointer value
 * NOTE: g_ctx must be set to call this function.  A global is required because
 * qsort doesn't take an opaque pointer.
 */
static int COVER_strict_cmp(const void *lp, const void *rp) {
  int result = COVER_cmp(g_ctx, lp, rp);
  if (result == 0) {
    result = lp < rp ? -1 : 1;
  }
  return result;
}

/**
 * Returns the first pointer in [first, last) whose element does not compare
 * less than value.  If no such element exists it returns last.
 */
static const size_t *COVER_lower_bound(const size_t *first, const size_t *last,
                                       size_t value) {
  size_t count = last - first;
  while (count != 0) {
    size_t step = count / 2;
    const size_t *ptr = first;
    ptr += step;
    if (*ptr < value) {
      first = ++ptr;
      count -= step + 1;
    } else {
      count = step;
    }
  }
  return first;
}

/**
 * Generic groupBy function.
 * Groups an array sorted by cmp into groups with equivalent values.
 * Calls grp for each group.
 */
static void
COVER_groupBy(const void *data, size_t count, size_t size, COVER_ctx_t *ctx,
              int (*cmp)(COVER_ctx_t *, const void *, const void *),
              void (*grp)(COVER_ctx_t *, const void *, const void *)) {
  const BYTE *ptr = (const BYTE *)data;
  size_t num = 0;
  while (num < count) {
    const BYTE *grpEnd = ptr + size;
    ++num;
    while (num < count && cmp(ctx, ptr, grpEnd) == 0) {
      grpEnd += size;
      ++num;
    }
    grp(ctx, ptr, grpEnd);
    ptr = grpEnd;
  }
}

static void COVER_group(COVER_ctx_t *ctx, const void *group,
                        const void *groupEnd) {
  /* The group consists of all the positions with the same first d bytes. */
  const U32 *grpPtr = (const U32 *)group;
  const U32 *grpEnd = (const U32 *)groupEnd;
  /* The dmerId is how we will reference this dmer.
   * This allows us to map the whole dmer space to a much smaller space, the
   * size of the suffix array.
   */
  const U32 dmerId = grpPtr - ctx->suffix;
  /* Count the number of samples this dmer shows up in */
  U32 freq = 0;
  /* Details */
  const size_t *curOffsetPtr = ctx->offsets;
  const size_t *offsetsEnd = ctx->offsets + ctx->nbSamples;
  /* Once *grpPtr >= curSampleEnd this occurrence of the dmer is in a
   * different sample than the last.
   */
  size_t curSampleEnd = ctx->offsets[0];
  for (; grpPtr != grpEnd; ++grpPtr) {
    /* Save the dmerId for this position so we can get back to it. */
    ctx->dmerAt[*grpPtr] = dmerId;
    /* Dictionaries only help for the first reference to the dmer.
     * After that zstd can reference the match from the previous reference.
     * So only count each dmer once for each sample it is in.
     */
    if (*grpPtr < curSampleEnd) {
      continue;
    }
    freq += 1;
    /* Binary search to find the end of the sample *grpPtr is in.
     * In the common case that grpPtr + 1 == grpEnd we can skip the binary
     * search because the loop is over.
     */
    if (grpPtr + 1 != grpEnd) {
      const size_t *sampleEndPtr =
          COVER_lower_bound(curOffsetPtr, offsetsEnd, *grpPtr);
      curSampleEnd = *sampleEndPtr;
      curOffsetPtr = sampleEndPtr + 1;
    }
  }
  /* At this point we are never going to look at this segment of the suffix
   * array again.  We take advantage of this fact to save memory.
   * We store the frequency of the dmer in the first position of the group,
   * which is dmerId.
   */
  ctx->suffix[dmerId] = freq;
}

/**
 * Selects the best segment in an epoch.
 * Segments of are scored according to the function:
 *
 * Let F(d) be the frequency of dmer d.
 * Let L(S) be the length of segment S.
 * Let S_i be the dmer at position i of segment S.
 *
 *                 F(S_1) + F(S_2) + ... + F(S_{L(S)-d+1})
 *     Score(S) = --------------------------------------
 *                          smoothing + L(S)
 *
 * We try each segment length in the range [kMin, kStep, kMax].
 * For each segment length we find the best segment according to Score.
 * We then take the best segment overall according to Score and return it.
 *
 * The difference from the paper is that we try multiple segment lengths.
 * We want to fit the segment length closer to the length of the useful part.
 * Longer segments allow longer matches, so they are worth more than shorter
 * ones.  However, if the extra length isn't high frequency it hurts us.
 * We add the smoothing in to give an advantage to longer segments.
 * The larger smoothing is, the more longer matches are favored.
 */
static COVER_segment_t COVER_selectSegment(COVER_ctx_t *ctx, U32 begin,
                                           U32 end) {
  const unsigned kMin = ctx->parameters.kMin;
  const unsigned kMax = ctx->parameters.kMax;
  const unsigned d = ctx->parameters.d;
  const unsigned smoothing = ctx->parameters.smoothing;
  /* Saves the best segment of any length tried */
  COVER_segment_t globalBestSegment = {0, 0, 0};
  /* For each segment length */
  U32 k;
  for (k = kMin; k <= kMax; k += ctx->parameters.kStep) {
    /* Save the best segment of this length */
    COVER_segment_t bestSegment = {0, 0, 0};
    COVER_segment_t activeSegment = {begin, begin, 0};
    const size_t dmersInK = k - d + 1;
    /* Reset the activeDmers in the segment */
    memset(ctx->activeDmers, 0, ctx->size * sizeof(*ctx->activeDmers));
    /* Slide the active segment through the whole epoch.
     * Save the best segment in bestSegment.
     */
    while (activeSegment.end < end) {
      /* The dmerId for the dmer at the next position */
      U32 newDmer = ctx->dmerAt[activeSegment.end];
      /* The entry in activeDmers for this dmerId */
      U32 *newDmerOcc = ctx->activeDmers + newDmer;
      /* If the dmer isn't already present in the segment add its score. */
      if (*newDmerOcc == 0) {
        /* The paper suggest using the L-0.5 norm, but experiments show that it
         * doesn't help.
         */
        activeSegment.score += ctx->freqs[newDmer];
      }
      /* Add the dmer to the segment */
      activeSegment.end += 1;
      *newDmerOcc += 1;

      /* If the window is now too large, drop the first position */
      if (activeSegment.end - activeSegment.begin == dmersInK + 1) {
        U32 delDmer = ctx->dmerAt[activeSegment.begin];
        U32 *delDmerOcc = ctx->activeDmers + delDmer;
        activeSegment.begin += 1;
        *delDmerOcc -= 1;
        /* If this is the last occurence of the dmer, subtract its score */
        if (*delDmerOcc == 0) {
          activeSegment.score -= ctx->freqs[delDmer];
        }
      }

      /* If this segment is the best so far save it */
      if (activeSegment.score > bestSegment.score) {
        bestSegment = activeSegment;
      }
    }
    {
      /* Trim off the zero frequency head and tail from the segment. */
      size_t newBegin = bestSegment.end;
      size_t newEnd = bestSegment.begin;
      size_t pos;
      for (pos = bestSegment.begin; pos != bestSegment.end; ++pos) {
        U32 freq = ctx->freqs[ctx->dmerAt[pos]];
        if (freq != 0) {
          newBegin = MIN(newBegin, pos);
          newEnd = pos + 1;
        }
      }
      bestSegment.begin = newBegin;
      bestSegment.end = newEnd;
      /* Calculate the final score normalizing for segment length */
      bestSegment.score /= (smoothing + (bestSegment.end - bestSegment.begin));
    }
    /* If this segment is the best so far for any length save it */
    if (bestSegment.score > globalBestSegment.score) {
      globalBestSegment = bestSegment;
    }
  }
  {
    /* Zero out the frequency of each dmer covered by the chosen segment. */
    size_t pos;
    for (pos = globalBestSegment.begin; pos != globalBestSegment.end; ++pos) {
      ctx->freqs[ctx->dmerAt[pos]] = 0;
    }
  }
  return globalBestSegment;
}

ZDICTLIB_API size_t COVER_trainFromBuffer(
    void *dictBuffer, size_t dictBufferCapacity, const void *samplesBuffer,
    const size_t *samplesSizes, size_t nbSamples, COVER_params_t parameters) {
  const size_t dictContentSize = dictBufferCapacity - COVER_ENTROPY_TABLES_SIZE;
  const size_t totalSamplesSize = COVER_sum(samplesSizes, nbSamples);
  BYTE *const dict = (BYTE *)dictBuffer;
  BYTE *const dictContent = dict + COVER_ENTROPY_TABLES_SIZE;
  const BYTE *const samples = (const BYTE *)samplesBuffer;
  /* Checks */
  if (nbSamples == 0) {
    return ERROR(GENERIC);
  }
  if (totalSamplesSize > (size_t)COVER_MAX_SAMPLES_SIZE) {
    return ERROR(GENERIC);
  }
  if (totalSamplesSize < parameters.d) {
    return ERROR(GENERIC);
  }
  if (parameters.d > parameters.kMin) {
    return ERROR(GENERIC);
  }
  if (parameters.kMin > parameters.kMax) {
    return ERROR(GENERIC);
  }
  if (dictBufferCapacity < COVER_ENTROPY_TABLES_SIZE) {
    return ERROR(dstSize_tooSmall);
  }
  g_displayLevel = parameters.notificationLevel;
  DISPLAYLEVEL(2, "Training on %zu samples of total size %zu\n", nbSamples,
               totalSamplesSize);
  {
    const size_t suffixSize = totalSamplesSize - parameters.d + 1;
    /* Partial suffix array */
    U32 *suffix = (U32 *)malloc(suffixSize * sizeof(U32));
    /* Maps index to the dmerID */
    U32 *dmerAt = (U32 *)malloc(suffixSize * sizeof(U32));
    /* Maps dmerID to the # of occurrences of the dmer in the segment */
    U32 *activeDmers = (U32 *)malloc(suffixSize * sizeof(U32));
    /* The offsets of each file */
    size_t *offsets = (size_t *)malloc((nbSamples + 1) * sizeof(size_t));
    size_t rc = 0;
    COVER_ctx_t ctx = {
        samples, offsets,     nbSamples,  suffix,     NULL,
        dmerAt,  activeDmers, suffixSize, parameters,
    };
    /* qsort doesn't take an opaque pointer, so pass the context as a global */
    g_ctx = &ctx;

    /* Fill offsets from the samlesSizes */
    {
      U32 i;
      offsets[0] = 0;
      for (i = 1; i <= nbSamples; ++i) {
        offsets[i] = offsets[i - 1] + samplesSizes[i - 1];
      }
    }
    DISPLAYLEVEL(2, "Constructing partial suffix array\n");
    {
      /* suffix is a partial suffix array.
       * It only sorts suffixes by their first parameters.d bytes.
       * The sort is stable, so each dmer group is sorted by position in input.
       */
      U32 i;
      for (i = 0; i < suffixSize; ++i) {
        suffix[i] = i;
      }
      qsort(suffix, suffixSize, sizeof(U32), &COVER_strict_cmp);
    }
    DISPLAYLEVEL(2, "Computing frequencies\n");
    /* For each dmer group (group of positions with the same first d bytes):
     * 1. For each position we set dmerAt[position] = dmerID.  The dmerID is
     *    (groupBeginPtr - suffix).  This allows us to go from position to
     *    dmerID so we can look up values in freq and activeDmers.
     * 2. We calculate how many samples the dmer occurs in and save it in
     *    freqs[dmerId].
     */
    COVER_groupBy(suffix, suffixSize, sizeof(U32), &ctx, &COVER_cmp,
                  &COVER_group);
    ctx.freqs = suffix;
    ctx.suffix = NULL;
    DISPLAYLEVEL(2, "Building dictionary\n");
    /* Select segments */
    {
      size_t tail = dictContentSize;
      /* Divide the data up into epochs of equal size.
       * We will select at least one segment from each epoch.
       */
      const U32 epochs = dictContentSize / parameters.kMax;
      const U32 epochSize = suffixSize / epochs;
      size_t epoch;
      DISPLAYLEVEL(3, "Breaking content into %u epochs of size %u\n", epochs,
                   epochSize);
      for (epoch = 0; tail > 0; epoch = (epoch + 1) % epochs) {
        const U32 epochBegin = epoch * epochSize;
        const U32 epochEnd = (epoch + 1) * epochSize;
        COVER_segment_t segment =
            COVER_selectSegment(&ctx, epochBegin, epochEnd);
        const size_t segmentSize =
            MIN(segment.end - segment.begin + parameters.d - 1, tail);
        if (segmentSize == 0) {
          break;
        }
        /* We fill the dictionary from the back to allow the best segments to be
         * referenced with the smallest offsets.
         */
        tail -= segmentSize;
        memcpy(dictContent + tail, samples + segment.begin, segmentSize);
        DISPLAYUPDATE(2, "\r%zu%%       ",
                      ((dictContentSize - tail) * 100) / dictContentSize);
      }
      {
        ZDICT_params_t zdictParams;
        memset(&zdictParams, 0, sizeof(zdictParams));
        zdictParams.notificationLevel = parameters.notificationLevel;
        zdictParams.dictID = parameters.dictID;
        zdictParams.compressionLevel = parameters.compressionLevel;
        rc = ZDICT_addEntropyTablesFromBuffer_advanced(
            dict, dictContentSize - tail, dictBufferCapacity, samplesBuffer,
            samplesSizes, nbSamples, zdictParams);
      }
    }
    DISPLAYLEVEL(2, "\r%79s\r", "");
    DISPLAYLEVEL(2, "Constructed dictionary of size %zu\n", rc);
    free(suffix);
    free(dmerAt);
    free(activeDmers);
    free(offsets);
    return rc;
  }
}
