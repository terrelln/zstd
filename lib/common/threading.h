
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
#ifndef WIN32_LEAN_AND_MEAN
#  define WIN32_LEAN_AND_MEAN
#endif

#include <windows.h>

/* mutex */
#define pthread_mutex_t           CRITICAL_SECTION
#define pthread_mutex_init(a,b)   InitializeCriticalSection((a))
#define pthread_mutex_destroy(a)  DeleteCriticalSection((a))
#define pthread_mutex_lock(a)     EnterCriticalSection((a))
#define pthread_mutex_unlock(a)   LeaveCriticalSection((a))

/**
 * Unless `ZSTD_CUSTOM_CV` or `ZSTD_NATIVE_CV` is already defined, detect Windows version.
 */
#if !defined(ZSTD_CUSTOM_CV) && !defined(ZSTD_NATIVE_CV)
#  if defined(WINVER) && WINVER >= 0x0600 && defined(_WIN32_WINNT) && _WIN32_WINNT >= 0x600
#    define ZSTD_NATIVE_CV
#  else
#    define ZSTD_CUSTOM_CV
#  endif
#endif

#if defined(ZSTD_NATIVE_CV)
/* Native condition variable for Windows Vista and newer */
#define pthread_cond_t             CONDITION_VARIABLE
#define pthread_cond_init(a, b)    InitializeConditionVariable((a))
#define pthread_cond_destroy(a)    /* No delete */
#define pthread_cond_wait(a, b)    SleepConditionVariableCS((a), (b), INFINITE)
#define pthread_cond_signal(a)     WakeConditionVariable((a))
#define pthread_cond_broadcast(a)  WakeAllConditionVariable((a))

#elif defined(ZSTD_CUSTOM_CV)
/* Custom condition variable for Windows XP and older */
#define pthread_cond_t             _pthread_cond_t
#define pthread_cond_init(a, b)    _pthread_cond_init((a), (b))
#define pthread_cond_destroy(a)    _pthread_cond_destroy((a))
#define pthread_cond_wait(a, b)    _pthread_cond_wait((a), (b))
#define pthread_cond_signal(a)     _pthread_cond_signal((a))
#define pthread_cond_broadcast(a)  _pthread_cond_broadcast((a))

typedef struct {
  pthread_mutex_t internalMutex;
  enum { SHUTDOWN = 0, RUNNING = 64212 } state;
  struct THREAD_event_s* waitingEvents;
  struct THREAD_event_s* recycledEvents;
} _pthread_cond_t;

int _pthread_cond_init(pthread_cond_t* cond, const void* attrs);
int _pthread_cond_destroy(pthread_cond_t* cond);
int _pthread_cond_wait(pthread_cond_t* cond, pthread_mutex_t* externalMutex);
int _pthread_cond_signal(pthread_cond_t* cond);
int _pthread_cond_broadcast(pthread_cond_t* cond);

#endif

/* pthread_create() and pthread_join() */
typedef struct {
    HANDLE handle;
    void* (*start_routine)(void*);
    void* arg;
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

#define pthread_mutex_t int   /* #define rather than typedef, as sometimes pthread support is implicit, resulting in duplicated symbols */
#define pthread_mutex_init(a,b)
#define pthread_mutex_destroy(a)
#define pthread_mutex_lock(a)
#define pthread_mutex_unlock(a)

#define pthread_cond_t int
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
