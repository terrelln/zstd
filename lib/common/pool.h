/**
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */
#ifndef POOL_H
#define POOL_H

#include <stddef.h>   /* size_t */

typedef struct POOL_ctx_s POOL_ctx;

/*! POOL_create() :
    Create a thread pool with at most `numThreads` threads.
    The maximum number of queued jobs before blocking is `2^queueLog - 1`.
    @return : The POOL_ctx on success else NULL.
*/
POOL_ctx *POOL_create(size_t numThreads, size_t queueLog);

/*! POOL_free() :
    Free a thread pool returned by POOL_create().
    @return : 0 on success and non-zero on failure.
*/
int POOL_free(POOL_ctx *ctx);

/*! POOL_function :
    The function type that can be added to a thread pool.
*/
typedef void (*POOL_function)(void *opaque);
/*! POOL_add_function :
    The function type for a generic thread pool add function.
    @returns : 0 on success and non-zero on error.
*/
typedef int (*POOL_add_function)(void *ctx, POOL_function function,
                                  void *opaque);

/*! POOL_add() :
    Add the job `function(opaque)` to the thread pool.
    Possibly blocks until there is room in the queue.
    @returns : 0 on success and non-zero on error.
*/
int POOL_add(POOL_ctx *ctx, POOL_function function, void *opaque);

/*! POOL_add_generic() :
    The same as POOL_add() except that it takes a void pointer for the context.
*/
int POOL_add_generic(void *ctx, POOL_function function, void *opaque);

#endif
