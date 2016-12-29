/**
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

/*-*************************************
*  Dependencies
***************************************/
#include <stdio.h>  /* fprintf */
#include <stdlib.h> /* malloc, free, qsort */
#include <string.h> /* memset */
#include <time.h>   /* clock */
#ifdef ZSTD_PTHREAD
#include <pthread.h>
#endif

#include "mem.h"           /* read */
#include "zstd_internal.h" /* includes zstd.h */
#ifndef ZDICT_STATIC_LINKING_ONLY
#define ZDICT_STATIC_LINKING_ONLY
#endif
#include "pool.h"
#include "zdict.h"
#include "zstd.h"

/*-*************************************
*  Constants
***************************************/
#define COVER_MAX_SAMPLES_SIZE ((U32)-1)

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
  size_t i;
  for (i = 0; i < nbSamples; ++i) {
    sum += samplesSizes[i];
  }
  return sum;
}

/**
 * A small specialized hash map for storing activeDmers.
 * The map does not resize, so if it becomes full it will loop forever.
 * Thus, the map must be large enough to store every value.
 * The map implements linear probing and keeps its load less than 0.5.
 */
#define MAP_EMPTY_VALUE ((U32)-1)
typedef struct COVER_map_pair_t_s {
  U32 key;
  U32 value;
} COVER_map_pair_t;

typedef struct COVER_map_s {
  COVER_map_pair_t *data;
  U32 sizeLog;
  U32 size;
  U32 sizeMask;
} COVER_map_t;

/**
 * Clear the map.
 */
static void COVER_map_clear(COVER_map_t *map) {
  memset(map->data, MAP_EMPTY_VALUE, map->size * sizeof(COVER_map_pair_t));
}

/**
 * Initializes a map of the given size.
 * Returns 1 on success and 0 on failure.
 * The map must be destroyed with COVER_map_destroy().
 * The map is only guaranteed to be large enough to hold size elements.
 */
static int COVER_map_init(COVER_map_t *map, U32 size) {
  map->sizeLog = ZSTD_highbit32(size - 1) + 2;
  map->size = (U32)1 << map->sizeLog;
  map->sizeMask = map->size - 1;
  map->data = (COVER_map_pair_t *)malloc(map->size * sizeof(COVER_map_pair_t));
  if (!map->data) {
    map->sizeLog = 0;
    map->size = 0;
    return 0;
  }
  COVER_map_clear(map);
  return 1;
}

/**
 * Internal hash function
 */
static const U32 prime4bytes = 2654435761U;
static U32 COVER_map_hash(COVER_map_t *map, U32 key) {
  return (key * prime4bytes) >> (32 - map->sizeLog);
}

/**
 * Helper function that returns the index that a key should be placed into.
 */
static U32 COVER_map_index(COVER_map_t *map, U32 key) {
  const U32 hash = COVER_map_hash(map, key);
  U32 i;
  for (i = hash;; i = (i + 1) & map->sizeMask) {
    COVER_map_pair_t *pos = &map->data[i];
    if (pos->value == MAP_EMPTY_VALUE) {
      return i;
    }
    if (pos->key == key) {
      return i;
    }
  }
}

/**
 * Returns the pointer to the value for key.
 * If key is not in the map, it is inserted and the value is set to 0.
 * The map must not be full.
 */
static U32 *COVER_map_at(COVER_map_t *map, U32 key) {
  COVER_map_pair_t *pos = &map->data[COVER_map_index(map, key)];
  if (pos->value == MAP_EMPTY_VALUE) {
    pos->key = key;
    pos->value = 0;
  }
  return &pos->value;
}

/**
 * Deletes key from the map if present.
 */
static void COVER_map_remove(COVER_map_t *map, U32 key) {
  U32 i = COVER_map_index(map, key);
  COVER_map_pair_t *del = &map->data[i];
  U32 shift = 1;
  if (del->value == MAP_EMPTY_VALUE) {
    return;
  }
  for (i = (i + 1) & map->sizeMask;; i = (i + 1) & map->sizeMask) {
    COVER_map_pair_t *const pos = &map->data[i];
    /* If the position is empty we are done */
    if (pos->value == MAP_EMPTY_VALUE) {
      del->value = MAP_EMPTY_VALUE;
      return;
    }
    /* If pos can be moved to del do so */
    if (((i - COVER_map_hash(map, pos->key)) & map->sizeMask) >= shift) {
      del->key = pos->key;
      del->value = pos->value;
      del = pos;
      shift = 1;
    } else {
      ++shift;
    }
  }
}

/**
 * Destroyes a map that is inited with COVER_map_init().
 */
static void COVER_map_destroy(COVER_map_t *map) {
  if (map->data) {
    free(map->data);
  }
  map->data = NULL;
  map->size = 0;
}

typedef struct {
  U32 begin;
  U32 end;
  double score;
} COVER_segment_t;

typedef struct {
  const BYTE *samples;
  size_t *offsets;
  const size_t *samplesSizes;
  size_t nbSamples;
  U32 *suffix;
  size_t suffixSize;
  U32 *freqs;
  U32 *dmerAt;
  unsigned d;
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
  return memcmp(ctx->samples + lhs, ctx->samples + rhs, ctx->d);
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
  const U32 dmerId = (U32)(grpPtr - ctx->suffix);
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
static COVER_segment_t COVER_selectSegment(const COVER_ctx_t *ctx, U32 *freqs,
                                           COVER_map_t *activeDmers, U32 begin,
                                           U32 end, COVER_params_t parameters) {
  /* Saves the best segment of any length tried */
  COVER_segment_t globalBestSegment = {0, 0, 0};
  /* For each segment length */
  U32 k;
  for (k = parameters.kMin; k <= parameters.kMax; k += parameters.kStep) {
    /* Save the best segment of this length */
    COVER_segment_t bestSegment = {0, 0, 0};
    COVER_segment_t activeSegment;
    const size_t dmersInK = k - ctx->d + 1;
    /* Reset the activeDmers in the segment */
    COVER_map_clear(activeDmers);
    activeSegment.begin = begin;
    activeSegment.end = begin;
    activeSegment.score = 0;
    /* Slide the active segment through the whole epoch.
     * Save the best segment in bestSegment.
     */
    while (activeSegment.end < end) {
      /* The dmerId for the dmer at the next position */
      U32 newDmer = ctx->dmerAt[activeSegment.end];
      /* The entry in activeDmers for this dmerId */
      U32 *newDmerOcc = COVER_map_at(activeDmers, newDmer);
      /* If the dmer isn't already present in the segment add its score. */
      if (*newDmerOcc == 0) {
        /* The paper suggest using the L-0.5 norm, but experiments show that it
         * doesn't help.
         */
        activeSegment.score += freqs[newDmer];
      }
      /* Add the dmer to the segment */
      activeSegment.end += 1;
      *newDmerOcc += 1;

      /* If the window is now too large, drop the first position */
      if (activeSegment.end - activeSegment.begin == dmersInK + 1) {
        U32 delDmer = ctx->dmerAt[activeSegment.begin];
        U32 *delDmerOcc = COVER_map_at(activeDmers, delDmer);
        activeSegment.begin += 1;
        *delDmerOcc -= 1;
        /* If this is the last occurence of the dmer, subtract its score */
        if (*delDmerOcc == 0) {
          COVER_map_remove(activeDmers, delDmer);
          activeSegment.score -= freqs[delDmer];
        }
      }

      /* If this segment is the best so far save it */
      if (activeSegment.score > bestSegment.score) {
        bestSegment = activeSegment;
      }
    }
    {
      /* Trim off the zero frequency head and tail from the segment. */
      U32 newBegin = bestSegment.end;
      U32 newEnd = bestSegment.begin;
      U32 pos;
      for (pos = bestSegment.begin; pos != bestSegment.end; ++pos) {
        U32 freq = freqs[ctx->dmerAt[pos]];
        if (freq != 0) {
          newBegin = MIN(newBegin, pos);
          newEnd = pos + 1;
        }
      }
      bestSegment.begin = newBegin;
      bestSegment.end = newEnd;
      /* Calculate the final score normalizing for segment length */
      bestSegment.score /=
          (parameters.smoothing + (bestSegment.end - bestSegment.begin));
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
      freqs[ctx->dmerAt[pos]] = 0;
    }
  }
  return globalBestSegment;
}

/**
 * Check the validity of the parameters.
 * If the parameters are valid and any are default, set them to the correct
 * values.
 * Returns 1 on success, 0 on failure.
 */
static int COVER_defaultParameters(COVER_params_t *parameters) {
  /* kMin and d are required parameters */
  if (parameters->d == 0 || parameters->kMin == 0) {
    return 0;
  }
  /* d <= kMin */
  if (parameters->d > parameters->kMin) {
    return 0;
  }
  /* If kMax is set (non-zero) then kMin <= kMax */
  if (parameters->kMax != 0 && parameters->kMax < parameters->kMin) {
    return 0;
  }
  /* If kMax is set, then kStep must be as well */
  if (parameters->kMax != 0 && parameters->kStep == 0) {
    return 0;
  }
  parameters->kMax = MAX(parameters->kMin, parameters->kMax);
  parameters->kStep = MAX(1, parameters->kStep);
  return 1;
}

static void COVER_destroy_ctx(COVER_ctx_t *ctx) {
  if (!ctx) {
    return;
  }
  if (ctx->suffix) {
    free(ctx->suffix);
  }
  if (ctx->freqs) {
    free(ctx->freqs);
  }
  if (ctx->dmerAt) {
    free(ctx->dmerAt);
  }
  if (ctx->offsets) {
    free(ctx->offsets);
  }
}

/**
 * Prepare a context for dictionary building.
 * The context is only dependent on the parameter `d` and can used multiple
 * times.
 * Returns 0 on success or non-zero on error.
 * The context must be destroyed with `COVER_destroy_ctx()`.
 */
static int COVER_init_ctx(COVER_ctx_t *ctx, const void *samplesBuffer,
                          const size_t *samplesSizes, unsigned nbSamples,
                          unsigned d) {
  const BYTE *const samples = (const BYTE *)samplesBuffer;
  const size_t totalSamplesSize = COVER_sum(samplesSizes, nbSamples);
  /* Checks */
  if (totalSamplesSize > (size_t)COVER_MAX_SAMPLES_SIZE) {
    return 1;
  }
  /* Zero the context */
  memset(ctx, 0, sizeof(*ctx));
  DISPLAYLEVEL(2, "Training on %u samples of total size %u\n", nbSamples,
               (U32)totalSamplesSize);
  ctx->samples = samples;
  ctx->samplesSizes = samplesSizes;
  ctx->nbSamples = nbSamples;
  /* Partial suffix array */
  ctx->suffixSize = totalSamplesSize - d + 1;
  ctx->suffix = (U32 *)malloc(ctx->suffixSize * sizeof(U32));
  /* Maps index to the dmerID */
  ctx->dmerAt = (U32 *)malloc(ctx->suffixSize * sizeof(U32));
  /* The offsets of each file */
  ctx->offsets = (size_t *)malloc((nbSamples + 1) * sizeof(size_t));
  if (!ctx->suffix || !ctx->dmerAt || !ctx->offsets) {
    COVER_destroy_ctx(ctx);
    return 1;
  }
  ctx->freqs = NULL;
  ctx->d = d;

  /* Fill offsets from the samlesSizes */
  {
    U32 i;
    ctx->offsets[0] = 0;
    for (i = 1; i <= nbSamples; ++i) {
      ctx->offsets[i] = ctx->offsets[i - 1] + samplesSizes[i - 1];
    }
  }
  DISPLAYLEVEL(2, "Constructing partial suffix array\n");
  {
    /* suffix is a partial suffix array.
     * It only sorts suffixes by their first parameters.d bytes.
     * The sort is stable, so each dmer group is sorted by position in input.
     */
    U32 i;
    for (i = 0; i < ctx->suffixSize; ++i) {
      ctx->suffix[i] = i;
    }
    /* qsort doesn't take an opaque pointer, so pass as a global */
    g_ctx = ctx;
    qsort(ctx->suffix, ctx->suffixSize, sizeof(U32), &COVER_strict_cmp);
  }
  DISPLAYLEVEL(2, "Computing frequencies\n");
  /* For each dmer group (group of positions with the same first d bytes):
   * 1. For each position we set dmerAt[position] = dmerID.  The dmerID is
   *    (groupBeginPtr - suffix).  This allows us to go from position to
   *    dmerID so we can look up values in freq.
   * 2. We calculate how many samples the dmer occurs in and save it in
   *    freqs[dmerId].
   */
  COVER_groupBy(ctx->suffix, ctx->suffixSize, sizeof(U32), ctx, &COVER_cmp,
                &COVER_group);
  ctx->freqs = ctx->suffix;
  ctx->suffix = NULL;
  return 0;
}

/**
 * Given the prepared context build the dictionary.
 */
static size_t COVER_buildDictionary(const COVER_ctx_t *ctx, U32 *freqs,
                                    COVER_map_t *activeDmers, void *dictBuffer,
                                    size_t dictBufferCapacity,
                                    COVER_params_t parameters) {
  BYTE *const dict = (BYTE *)dictBuffer;
  size_t tail = dictBufferCapacity;
  /* Divide the data up into epochs of equal size.
   * We will select at least one segment from each epoch.
   */
  const U32 epochs = (U32)(dictBufferCapacity / parameters.kMax);
  const U32 epochSize = (U32)(ctx->suffixSize / epochs);
  size_t epoch;
  DISPLAYLEVEL(3, "Breaking content into %u epochs of size %u\n", epochs,
               epochSize);
  for (epoch = 0; tail > 0; epoch = (epoch + 1) % epochs) {
    const U32 epochBegin = (U32)(epoch * epochSize);
    const U32 epochEnd = epochBegin + epochSize;
    size_t segmentSize;
    COVER_segment_t segment = COVER_selectSegment(
        ctx, freqs, activeDmers, epochBegin, epochEnd, parameters);
    segmentSize = MIN(segment.end - segment.begin + parameters.d - 1, tail);
    if (segmentSize == 0) {
      break;
    }
    /* We fill the dictionary from the back to allow the best segments to be
     * referenced with the smallest offsets.
     */
    tail -= segmentSize;
    memcpy(dict + tail, ctx->samples + segment.begin, segmentSize);
    DISPLAYUPDATE(
        3, "\r%u%%       ",
        (U32)(((dictBufferCapacity - tail) * 100) / dictBufferCapacity));
  }
  DISPLAYLEVEL(3, "\r%79s\r", "");
  return tail;
}

/**
 * Translate from COVER_params_t to ZDICT_params_t required for finalizing the
 * dictionary.
 */
static ZDICT_params_t COVER_translateParams(COVER_params_t parameters) {
  ZDICT_params_t zdictParams;
  memset(&zdictParams, 0, sizeof(zdictParams));
  zdictParams.notificationLevel = parameters.notificationLevel;
  zdictParams.dictID = parameters.dictID;
  zdictParams.compressionLevel = parameters.compressionLevel;
  return zdictParams;
}

/**
 * Constructs a dictionary using a heuristic based on the following paper:
 *
 * Liao, Petri, Moffat, Wirth
 * Effective Construction of Relative Lempel-Ziv Dictionaries
 * Published in WWW 2016.
 */
ZDICTLIB_API size_t COVER_trainFromBuffer(
    void *dictBuffer, size_t dictBufferCapacity, const void *samplesBuffer,
    const size_t *samplesSizes, unsigned nbSamples, COVER_params_t parameters) {
  BYTE *const dict = (BYTE *)dictBuffer;
  COVER_ctx_t ctx;
  COVER_map_t activeDmers;
  size_t rc;
  /* Checks */
  if (!COVER_defaultParameters(&parameters)) {
    DISPLAYLEVEL(1, "Cover parameters incorrect\n");
    return ERROR(GENERIC);
  }
  if (nbSamples == 0) {
    DISPLAYLEVEL(1, "Cover must have at least one input file\n");
    return ERROR(GENERIC);
  }
  if (dictBufferCapacity < ZDICT_DICTSIZE_MIN) {
    return ERROR(dstSize_tooSmall);
  }
  /* Initialize global data */
  g_displayLevel = parameters.notificationLevel;
  /* Initialize context and activeDmers */
  if (COVER_init_ctx(&ctx, samplesBuffer, samplesSizes, nbSamples,
                     parameters.d)) {
    DISPLAYLEVEL(1, "Failed to initialize context\n");
    rc = ERROR(GENERIC);
    goto _cleanup;
  }
  if (!COVER_map_init(&activeDmers, parameters.kMax - parameters.d + 1)) {
    DISPLAYLEVEL(1, "Failed to allocate dmer map: out of memory\n");
    rc = ERROR(GENERIC);
    goto _cleanup;
  }

  DISPLAYLEVEL(2, "Building dictionary\n");
  {
    const size_t tail =
        COVER_buildDictionary(&ctx, ctx.freqs, &activeDmers, dictBuffer,
                              dictBufferCapacity, parameters);
    ZDICT_params_t zdictParams = COVER_translateParams(parameters);
    DISPLAYLEVEL(2, "Dictionary content size: %u",
                 (U32)(dictBufferCapacity - tail));
    rc = ZDICT_finalizeDictionary(dict, dictBufferCapacity, dict + tail,
                                  dictBufferCapacity - tail, samplesBuffer,
                                  samplesSizes, nbSamples, zdictParams);
  }
  if (!ZSTD_isError(rc)) {
    DISPLAYLEVEL(2, "Constructed dictionary of size %u\n", (U32)rc);
  }
_cleanup:
  COVER_destroy_ctx(&ctx);
  COVER_map_destroy(&activeDmers);
  return rc;
}

typedef struct COVER_best_s {
#ifdef ZSTD_PTHREAD
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  size_t liveJobs;
#endif
  COVER_params_t parameters;
  size_t size;
} COVER_best_t;

static int COVER_best_init(COVER_best_t *best) {
  if (!best) {
    return 1;
  }
#ifdef ZSTD_PTHREAD
  if (pthread_mutex_init(&best->mutex, NULL)) {
    return 1;
  }
  if (pthread_cond_init(&best->cond, NULL)) {
    return 1;
  }
  best->liveJobs = 0;
#endif
  best->size = (size_t)-1;
  memset(&best->parameters, 0, sizeof(best->parameters));
  return 0;
}

static int COVER_best_wait(COVER_best_t *best) {
  int rc = 0;
  if (!best) {
    return 1;
  }
#ifdef ZSTD_PTHREAD
  rc |= pthread_mutex_lock(&best->mutex);
  while (best->liveJobs != 0 && rc == 0) {
    rc |= pthread_cond_wait(&best->cond, &best->mutex);
  }
  rc |= pthread_mutex_unlock(&best->mutex);
#endif
  return rc;
}

static int COVER_best_destroy(COVER_best_t *best) {
  int rc = 0;
  if (!best) {
    return 0;
  }
  rc |= COVER_best_wait(best);
#ifdef ZSTD_PTHREAD
  rc |= pthread_mutex_destroy(&best->mutex);
  rc |= pthread_cond_destroy(&best->cond);
#endif
  return rc;
}

static int COVER_best_start(COVER_best_t *best) {
  int rc = 0;
  if (!best) {
    return 1;
  }
#ifdef ZSTD_PTHREAD
  rc |= pthread_mutex_lock(&best->mutex);
  ++best->liveJobs;
  rc |= pthread_mutex_unlock(&best->mutex);
#endif
  return rc;
}

static int COVER_best_finish(COVER_best_t *best, COVER_params_t parameters,
                             size_t size) {
  int rc = 0;
  if (!best) {
    return 1;
  }
  {
#ifdef ZSTD_PTHREAD
    size_t liveJobs = best->liveJobs;
    rc |= pthread_mutex_lock(&best->mutex);
    --best->liveJobs;
#endif
    if (size < best->size) {
      best->parameters = parameters;
      best->size = size;
    }
#ifdef ZSTD_PTHREAD
    rc |= pthread_mutex_unlock(&best->mutex);
    if (liveJobs == 0) {
      pthread_cond_broadcast(&best->cond);
    }
#endif
  }
  return rc;
}

typedef struct COVER_tryParameters_data_s {
  const COVER_ctx_t *ctx;
  COVER_best_t *best;
  size_t dictBufferCapacity;
  COVER_params_t parameters;
} COVER_tryParameters_data_t;

static void COVER_tryParameters(void *opaque) {
  COVER_tryParameters_data_t *data = (COVER_tryParameters_data_t *)opaque;
  const COVER_ctx_t *ctx = data->ctx;
  COVER_params_t parameters = data->parameters;
  size_t dictBufferCapacity = data->dictBufferCapacity;
  COVER_map_t activeDmers;
  BYTE *const dict = (BYTE * const)malloc(dictBufferCapacity);
  U32 *freqs = (U32 *)malloc(ctx->suffixSize * sizeof(U32));
  if (!COVER_map_init(&activeDmers, parameters.kMax - parameters.d + 1)) {
    DISPLAYLEVEL(1, "Failed to allocate dmer map: out of memory\n");
    goto _cleanup;
  }
  if (!dict || !freqs) {
    DISPLAYLEVEL(1, "Failed to allocate dictionary buffer\n");
    goto _cleanup;
  }
  memcpy(freqs, ctx->freqs, ctx->suffixSize * sizeof(U32));
  {
    const size_t tail = COVER_buildDictionary(ctx, freqs, &activeDmers, dict,
                                              dictBufferCapacity, parameters);
    ZDICT_params_t zdictParams = COVER_translateParams(parameters);
    dictBufferCapacity = ZDICT_finalizeDictionary(
        dict, dictBufferCapacity, dict + tail, dictBufferCapacity - tail,
        ctx->samples, ctx->samplesSizes, ctx->nbSamples, zdictParams);
    if (ZDICT_isError(dictBufferCapacity)) {
      DISPLAYLEVEL(1, "Failed to finalize dictionary\n");
      goto _cleanup;
    }
  }
  /* Check total compressed size */
  {
    /* Pointers */
    ZSTD_CCtx *cctx;
    ZSTD_CDict *cdict;
    void *dst;
    /* Local variables */
    size_t totalCompressedSize = 0;
    size_t dstCapacity;
    size_t i;
    /* Allocate dst with enough space to compress the maximum sized sample */
    {
      size_t maxSampleSize = 0;
      for (i = 0; i < ctx->nbSamples; ++i) {
        maxSampleSize = MAX(ctx->samplesSizes[i], maxSampleSize);
      }
      dstCapacity = ZSTD_compressBound(maxSampleSize);
      dst = malloc(dstCapacity);
    }
    /* Create the cctx and cdict */
    cctx = ZSTD_createCCtx();
    cdict =
        ZSTD_createCDict(dict, dictBufferCapacity, parameters.compressionLevel);
    if (!dst || !cctx || !cdict) {
      goto _compressError;
    }
    /* Compress each sample and sum their sizes (or error) */
    for (i = 0; i < ctx->nbSamples; ++i) {
      const size_t size = ZSTD_compress_usingCDict(
          cctx, dst, dstCapacity, ctx->samples + ctx->offsets[i],
          ctx->samplesSizes[i], cdict);
      if (ZSTD_isError(size)) {
        goto _compressError;
      }
      totalCompressedSize += size;
    }
    COVER_best_finish(data->best, parameters, totalCompressedSize);
    goto _cleanup;
  _compressError:
    COVER_best_finish(data->best, parameters, ERROR(GENERIC));
    ZSTD_freeCCtx(cctx);
    ZSTD_freeCDict(cdict);
    if (dst) {
      free(dst);
    }
  }

_cleanup:
  COVER_map_destroy(&activeDmers);
  free(data);
  if (dict) {
    free(dict);
  }
  if (freqs) {
    free(freqs);
  }
}

ZDICTLIB_API size_t COVER_optimizeParameters(
    size_t dictBufferCapacity, const void *samplesBuffer,
    const size_t *samplesSizes, unsigned nbSamples, void *threadPool,
    POOL_add_function addJob, COVER_params_t *parameters) {
  unsigned d;
  COVER_best_t best;
  if (COVER_best_init(&best)) {
    return ERROR(GENERIC);
  }
  g_displayLevel = parameters->notificationLevel;
  for (d = 6; d <= 12; d += 2) {
    COVER_ctx_t ctx;
    const unsigned min = 32;
    const unsigned max = 1024;
    unsigned kMin = min;
    unsigned kMax = max;
    DISPLAYLEVEL(2, "d=%u\n", d);
    if (COVER_init_ctx(&ctx, samplesBuffer, samplesSizes, nbSamples, d)) {
      DISPLAYLEVEL(1, "Failed to initialize context\n");
      COVER_best_destroy(&best);
      return ERROR(GENERIC);
    }
    for (kMin = min; kMin <= max; kMin <<= 1) {
      DISPLAYLEVEL(2, "kMin=%u\n", kMin);
      for (kMax = kMin; kMax <= max; kMax <<= 1) {
        unsigned kStep = MAX((kMax - kMin) >> 3, 1);
        unsigned smoothing;
        DISPLAYLEVEL(2, "kMax=%u\n", kMax);
        for (smoothing = kMin >> 2; smoothing <= (kMin << 1); smoothing <<= 1) {
          COVER_tryParameters_data_t *data =
              malloc(sizeof(COVER_tryParameters_data_t));
          DISPLAYLEVEL(2, "smoothing=%u\n", smoothing);
          if (!data) {
            DISPLAYLEVEL(1, "Failed to allocate parameters\n");
            COVER_best_destroy(&best);
            COVER_destroy_ctx(&ctx);
            return ERROR(GENERIC);
          }
          data->ctx = &ctx;
          data->best = &best;
          data->dictBufferCapacity = dictBufferCapacity;
          data->parameters = *parameters;
          data->parameters.d = d;
          data->parameters.kMin = kMin;
          data->parameters.kStep = kStep;
          data->parameters.kMax = kMax;
          data->parameters.smoothing = smoothing;
          COVER_best_start(&best);
          addJob(threadPool, &COVER_tryParameters, data);
        }
      }
    }
    COVER_best_wait(&best);
    COVER_destroy_ctx(&ctx);
  }
  {
    const size_t size = best.size;
    *parameters = best.parameters;
    COVER_best_destroy(&best);
    return size;
  }
}
