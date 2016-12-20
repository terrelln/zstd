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
typedef unsigned COVER_score_t;
/* The norm to use when adding frequencies */
#define COVER_NORM(frequency) (frequency)
/* Normalize the score across all segment lengths */
#define COVER_NORMALIZE_SCORE(segment)                                         \
  (COVER_score_t)((segment).score / pow((segment).end - (segment).begin, 0.9))

/*-*************************************
*  Dependencies
***************************************/
#include <math.h>   /* pow */
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

/*-*************************************
*  Console display
***************************************/
// #define display(...)                                                           \
//   {                                                                            \
//     fprintf(stderr, __va_args__);                                              \
//     fflush(stderr);                                                            \
//   }
// #define displaylevel(l, ...)                                                   \
//   if (notificationlevel >= l) {                                                \
//     display(__va_args__);                                                      \
//   } /* 0 : no display;   1: errors;   2: default;  3: details;  4: debug */
//
// static clock_t zdict_clockspan(clock_t nprevious) {
//   return clock() - nprevious;
// }

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
  COVER_params_t parameters;
} COVER_ctx_t;

static int COVER_cmp(void *opaque, const void *lp, const void *rp) {
  COVER_ctx_t *ctx = (COVER_ctx_t *)opaque;
  const U32 lhs = *(const U32 *)lp;
  const U32 rhs = *(const U32 *)rp;
  return memcmp(ctx->samples + lhs, ctx->samples + rhs, ctx->parameters.cover);
}

static int COVER_strict_cmp(void *opaque, const void *lp, const void *rp) {
  int result = COVER_cmp(opaque, lp, rp);
  if (result == 0) {
    result = lp < rp ? -1 : 1;
  }
  return result;
}

static const size_t *COVER_lower_bound(const size_t *first, const size_t *last,
                                       size_t value) {
  size_t count = last - first;
  size_t step;
  while (count != 0) {
    const size_t *ptr = first;
    step = count / 2;
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

static void COVER_group(void *opaque, const void *group, const void *groupEnd) {
  COVER_ctx_t *ctx = (COVER_ctx_t *)opaque;
  const U32 *grpPtr = (const U32 *)group;
  const U32 *grpEnd = (const U32 *)groupEnd;
  const U32 dmerId = grpPtr - ctx->suffix;
  /* Count the number of samples this dmer shows up in */
  U32 freq = 0;
  const size_t *curOffsetPtr = ctx->offsets;
  const size_t *offsetsEnd = ctx->offsets + ctx->nbSamples;
  size_t curSampleEnd = ctx->offsets[0];
  for (; grpPtr != grpEnd; ++grpPtr) {
    ctx->dmerAt[*grpPtr] = dmerId;
    if (*grpPtr < curSampleEnd) {
      continue;
    }
    freq += 1;
    if (grpPtr + 1 != grpEnd) {
      const size_t *sampleEndPtr =
          COVER_lower_bound(curOffsetPtr, offsetsEnd, *grpPtr + 1);
      curSampleEnd = *sampleEndPtr;
      curOffsetPtr = sampleEndPtr + 1;
    }
  }
  ctx->freqs[dmerId] = freq;
}

static void COVER_groupBy(const void *data, size_t count, size_t size,
                          void *opaque,
                          int (*cmp)(void *, const void *, const void *),
                          void (*grp)(void *, const void *, const void *)) {
  /* TODO: Parallelize */
  const BYTE *ptr = (const BYTE *)data;
  size_t num = 0;
  while (num < count) {
    const BYTE *grpEnd = ptr + size;
    ++num;
    while (num < count && cmp(opaque, ptr, grpEnd) == 0) {
      grpEnd += size;
      ++num;
    }
    grp(opaque, ptr, grpEnd);
    ptr = grpEnd;
  }
}

static COVER_segment_t COVER_selectSegment(COVER_ctx_t *ctx, U32 begin,
                                           U32 end) {
  const unsigned minSegment = ctx->parameters.minSegment;
  const unsigned maxSegment = ctx->parameters.maxSegment;
  const unsigned cover = ctx->parameters.cover;
  COVER_segment_t globalBestSegment = {0, 0, 0};
  U32 segmentSize;
  for (segmentSize = minSegment; segmentSize <= maxSegment;
       segmentSize += ctx->parameters.step) {
    COVER_segment_t bestSegment = {0, 0, 0};
    COVER_segment_t activeSegment = {begin, begin, 0};
    const size_t dmersInSegment = segmentSize - cover + 1;
    memset(ctx->activeDmers, 0, (U32)1 << cover);

    while (activeSegment.end < end) {
      U32 newDmer = ctx->dmerAt[activeSegment.end];
      U32 *newDmerOcc = ctx->activeDmers + newDmer;
      activeSegment.end += 1;
      if (*newDmerOcc == 0) {
        activeSegment.score += COVER_NORM(ctx->freqs[newDmer]);
      }
      *newDmerOcc += 1;

      if (activeSegment.end - activeSegment.begin == dmersInSegment + 1) {
        U32 delDmer = ctx->dmerAt[activeSegment.begin];
        U32 *delDmerOcc = ctx->activeDmers + delDmer;
        activeSegment.begin += 1;
        *delDmerOcc -= 1;
        if (*delDmerOcc == 0) {
          activeSegment.score -= COVER_NORM(ctx->freqs[delDmer]);
        }
      }

      if (activeSegment.score > bestSegment.score) {
        bestSegment = activeSegment;
      }
    }
    {
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
      bestSegment.score = COVER_NORMALIZE_SCORE(bestSegment);
    }
    if (bestSegment.score > globalBestSegment.score) {
      globalBestSegment = bestSegment;
    }
  }
  {
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
  const size_t totalSamplesSize = COVER_sum(samplesSizes, nbSamples);
  BYTE *const dict = (BYTE *)dictBuffer;
  const BYTE *const samples = (const BYTE *)samplesBuffer;
  /* Checks */
  if (nbSamples == 0) {
    return ERROR(GENERIC);
  }
  if (totalSamplesSize > (size_t)COVER_MAX_SAMPLES_SIZE) {
    return ERROR(GENERIC);
  }
  if (parameters.cover > parameters.minSegment) {
    return ERROR(GENERIC);
  }
  {
    const size_t suffixSize = totalSamplesSize - parameters.cover + 1;
    U32 *suffix = (U32 *)malloc(suffixSize * sizeof(U32));
    U32 *dmerAt = (U32 *)malloc(suffixSize * sizeof(U32));
    U32 *activeDmers = (U32 *)malloc(suffixSize * sizeof(U32));
    size_t *offsets = (size_t *)malloc(nbSamples * sizeof(size_t));
    size_t rc = 0;
    COVER_ctx_t ctx = {
        (const BYTE *)samplesBuffer,
        offsets,
        nbSamples,
        suffix,
        NULL,
        dmerAt,
        activeDmers,
        parameters,
    };

    /* Fill offsets */
    {
      U32 i;
      offsets[0] = 0;
      for (i = 1; i < nbSamples; ++i) {
        offsets[i] = samplesSizes[i - 1];
      }
    }
    /* Construct partial suffix array */
    {
      U32 i;
      for (i = 0; i < suffixSize; ++i) {
        suffix[i] = i;
      }
      qsort_r(suffix, suffixSize, sizeof(U32), &ctx, &COVER_strict_cmp);
    }
    /* Compute frequencies for each dmer */
    ctx.freqs = suffix;
    COVER_groupBy(suffix, suffixSize, sizeof(U32), &ctx, &COVER_cmp,
                  &COVER_group);
    ctx.suffix = NULL;
    /* Select segments */
    {
      size_t tail = dictBufferCapacity;
      const U32 epochs = dictBufferCapacity / parameters.maxSegment;
      const U32 epochSize = suffixSize / epochs;
      size_t epoch;
      for (epoch = 0; tail > 0; epoch = (epoch + 1) % epochs) {
        const U32 epochBegin = epoch * epochSize;
        const U32 epochEnd = (epoch + 1) * epochSize;
        COVER_segment_t segment =
            COVER_selectSegment(&ctx, epochBegin, epochEnd);
        const size_t segmentSize = MIN(segment.end - segment.begin, tail);
        if (segmentSize == 0) {
          break;
        }
        tail -= segmentSize;
        memcpy(dict + tail, samples + segment.begin, segmentSize);
      }
      if (tail != 0) {
        rc = dictBufferCapacity - tail;
        memmove(dict, dict + tail, rc);
      }
    }

    free(suffix);
    free(dmerAt);
    free(activeDmers);
    free(offsets);
    return rc;
  }
}
