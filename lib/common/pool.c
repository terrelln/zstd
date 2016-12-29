/**
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#include "pool.h"
#include <stddef.h>  /* size_t */
#include <stdlib.h>  /* malloc, calloc, free */

#ifdef ZSTD_PTHREAD

#include <pthread.h>
typedef struct POOL_job_s {
  POOL_function function;
  void *opaque;
} POOL_job;

struct POOL_ctx_s {
    pthread_t *threads;
    size_t numThreads;

    POOL_job *queue;
    size_t queueHead;
    size_t queueTail;
    size_t queueMask;
    pthread_mutex_t queueMutex;
    pthread_cond_t queuePushCond;
    pthread_cond_t queuePopCond;
    int shutdown;
};

static void *POOL_thread(void *opaque) {
    POOL_ctx *ctx = (POOL_ctx *)opaque;
    if (!ctx) { return NULL; }
    for (;;) {
        /* Lock the mutex and wait for a non-empty queue or until shutdown */
        if (pthread_mutex_lock(&ctx->queueMutex)) { return NULL; }
        while (ctx->queueHead == ctx->queueTail && !ctx->shutdown) {
            if (pthread_cond_wait(&ctx->queuePopCond, &ctx->queueMutex)) { return NULL; }
        }
        /* empty => shutting down: so stop */
        if (ctx->queueHead == ctx->queueTail) {
            if (pthread_mutex_unlock(&ctx->queueMutex)) { return NULL; }
            return opaque;
        }
        {
            /* Pop a job off the queue */
            POOL_job job = ctx->queue[ctx->queueHead];
            ctx->queueHead = (ctx->queueHead + 1) & ctx->queueMask;
            /* Unlock the mutex, signal a pusher, and run the job */
            if (pthread_mutex_unlock(&ctx->queueMutex)) { return NULL; }
            if (pthread_cond_signal(&ctx->queuePushCond)) { return NULL; }
            job.function(job.opaque);
        }
    }
    /* Unreachable */
}

POOL_ctx *POOL_create(size_t numThreads, size_t queueLog) {
    size_t const queueSize = (size_t)1 << queueLog;
    int err = 0;
    /* Allocate the context and zero initialize */
    POOL_ctx *ctx = (POOL_ctx *)calloc(1, sizeof(POOL_ctx));
    if (!ctx) { return NULL; }
    /* Initialize the job queue */
    ctx->queue = (POOL_job *)malloc(queueSize * sizeof(POOL_job));
    ctx->queueHead = 0;
    ctx->queueTail = 0;
    ctx->queueMask = queueSize - 1;
    err |= pthread_mutex_init(&ctx->queueMutex, NULL);
    err |= pthread_cond_init(&ctx->queuePushCond, NULL);
    err |= pthread_cond_init(&ctx->queuePopCond, NULL);
    ctx->shutdown = 0;
    /* Allocate space for the thread handles */
    ctx->threads = (pthread_t *)malloc(numThreads * sizeof(pthread_t));
    ctx->numThreads = 0;
    /* Check for errors */
    if (!ctx->threads || !ctx->queue || err) { POOL_free(ctx); return NULL; }
    /* Initialize the threads */
    {   size_t i;
        for (i = 0; i < numThreads; ++i) {
            if (pthread_create(&ctx->threads[i], NULL, &POOL_thread, ctx)) {
                ctx->numThreads = i;
                POOL_free(ctx);
                return NULL;
            }
        }
        ctx->numThreads = numThreads;
    }
    return ctx;
}

/*! POOL_join() :
    Join all of the threads.
    Returns non-zero if any of the joins failed, or if any thread returned an
    error.
    Returns zero on success.
*/
static int POOL_join(POOL_ctx *ctx) {
    int err = 0;
    /* Shut down the queue */
    err |= pthread_mutex_lock(&ctx->queueMutex);
    ctx->shutdown = 1;
    err |= pthread_mutex_unlock(&ctx->queueMutex);
    /* Wake up sleeping threads */
    err |= pthread_cond_broadcast(&ctx->queuePushCond);
    err |= pthread_cond_broadcast(&ctx->queuePopCond);
    /* Join all of the threads */
    {   size_t i;
        for (i = 0; i < ctx->numThreads; ++i) {
            void *retval;
            err |= pthread_join(ctx->threads[i], &retval);
            err |= !retval;
        }
    }
    return err;
}

int POOL_free(POOL_ctx *ctx) {
    int err = 0;
    if (!ctx) { return 0; }
    err = POOL_join(ctx);
    err |= pthread_mutex_destroy(&ctx->queueMutex);
    err |= pthread_cond_destroy(&ctx->queuePushCond);
    err |= pthread_cond_destroy(&ctx->queuePopCond);
    if (ctx->queue) free(ctx->queue);
    if (ctx->threads) free(ctx->threads);
    free(ctx);
    return err;
}

int POOL_add(POOL_ctx *ctx, POOL_function function, void *opaque) {
    int err = 0;
    if (!ctx) { return 1; }

    err |= pthread_mutex_lock(&ctx->queueMutex);
    {
        POOL_job job = {function, opaque};
        size_t newTail = (ctx->queueTail + 1) & ctx->queueMask;
        while (ctx->queueHead == newTail && !(ctx->shutdown || err)) {
          err |= pthread_cond_wait(&ctx->queuePushCond, &ctx->queueMutex);
          newTail = (ctx->queueTail + 1) & ctx->queueMask;
        }
        if (!(ctx->shutdown || err)) {
            ctx->queue[ctx->queueTail] = job;
            ctx->queueTail = newTail;
        }
    }
    err |= pthread_mutex_unlock(&ctx->queueMutex);
    err |= pthread_cond_signal(&ctx->queuePopCond);
    return err;
}

#else

struct POOL_ctx_s {
  int _;
};

POOL_ctx *POOL_create(size_t numThreads, size_t queueLog) {
  (void)numThreads;
  (void)queueLog;
  return (POOL_ctx *)malloc(sizeof(POOL_ctx));
}

int POOL_free(POOL_ctx *ctx) {
  if (ctx) free(ctx);
  return 0;
}

int POOL_add(POOL_ctx *ctx, POOL_function function, void *opaque) {
  (void)ctx;
  function(opaque);
  return 0;
}

#endif

int POOL_add_generic(void *ctx, POOL_function function, void *opaque) {
  return POOL_add((POOL_ctx *)ctx, function, opaque);
}
