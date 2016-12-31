
/**
 * Copyright (c) 2016 Tino Reichardt
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 *
 * You can contact the author at:
 * - zstdmt source repository: https://github.com/mcmilk/zstdmt
 */

#ifndef THREADING_H_938743
#define THREADING_H_938743

#if defined (__cplusplus)
extern "C" {
#endif

#if defined(ZSTD_PTHREAD) && defined(_WIN32)

/**
 * Windows minimalist Pthread Wrapper, based on :
 * http://www.cse.wustl.edu/~schmidt/win32-cv-1.html
 */

#include <windows.h>

/* mutex */
#define pthread_mutex_t           CRITICAL_SECTION
#define pthread_mutex_init(a,b)   InitializeCriticalSection((a))
#define pthread_mutex_destroy(a)  DeleteCriticalSection((a))
#define pthread_mutex_lock(a)     EnterCriticalSection((a))
#define pthread_mutex_unlock(a)   LeaveCriticalSection((a))

/* condition variable */
#define pthread_cond_t             CONDITION_VARIABLE
#define pthread_cond_init(a, b)    InitializeConditionVariable((a))
#define pthread_cond_destroy(a)    /* No delete */
#define pthread_cond_wait(a, b)    SleepConditionVariableCS((a), (b), INFINITE)
#define pthread_cond_signal(a)     WakeConditionVariable((a))
#define pthread_cond_broadcast(a)  WakeAllConditionVariable((a))

/* pthread_create() and pthread_join() */
typedef struct {
    HANDLE handle;
    void* (*start_routine)(void*);
    void*varg;
} pthread_t;

int pthread_create(pthread_t* thread, const void* unused,
                   void* (*start_routine) (void*), void* arg);

#define pthread_join(a, b) _pthread_join(&(a), (b))
int _pthread_join(pthread_t* thread, void** value_ptr);

/**
 * add here more wrappers as required
 */


#elif defined(ZSTD_PTHREAD)   /* posix assumed ; need a better detection mathod */
/* ===   POSIX Systems   === */
#  include <pthread.h>

#else  /* ZSTD_PTHREAD not defined */
/* No multithreading support */

typedef int pthread_mutex_t;
#define pthread_mutex_init(a,b)
#define pthread_mutex_destroy(a)
#define pthread_mutex_lock(a)
#define pthread_mutex_unlock(a)

typedef int pthread_cond_t;
#define pthread_cond_init(a,b)
#define pthread_cond_destroy(a)
#define pthread_cond_wait(a,b)
#define pthread_cond_signal(a)
#define pthread_cond_broadcast(a)

/* do not use pthread_t */

#endif /* ZSTD_PTHREAD */

#if defined (__cplusplus)
}
#endif

#endif /* THREADING_H_938743 */
