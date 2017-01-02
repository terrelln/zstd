
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

/**
 * This file will hold wrapper for systems, which do not support Pthreads
 */

/**
 * Compiler specifics
 */
#if defined(__GNUC__)
#  define DEBUG_ONLY __attribute__((unused))
#else
#  define DEBUG_ONLY
#endif

#if defined(ZSTD_PTHREAD) && defined(_WIN32)

/**
 * Windows minimalist Pthread Wrapper, based on :
 * http://www.cse.wustl.edu/~schmidt/win32-cv-1.html
 */


/* ===  Dependencies  === */
#include <process.h>
#include <errno.h>
#include "threading.h"


/* === Debugging === */
#ifndef NDEBUG
#  ifndef ZSTD_DEBUG
#    define NDEBUG
#  endif
#endif
#include <assert.h>


/* ===  Implementation  === */

static unsigned __stdcall worker(void *arg)
{
    pthread_t* const thread = (pthread_t*) arg;
    thread->arg = thread->start_routine(thread->arg);
    return 0;
}

int pthread_create(pthread_t* thread, const void* unused,
            void* (*start_routine) (void*), void* arg)
{
    (void)unused;
    thread->arg = arg;
    thread->start_routine = start_routine;
    thread->handle = (HANDLE) _beginthreadex(NULL, 0, worker, thread, 0, NULL);

    if (!thread->handle)
        return errno;
    else
        return 0;
}

int _pthread_join(pthread_t * thread, void **value_ptr)
{
    DWORD result;

    if (!thread->handle) return 0;

    result = WaitForSingleObject(thread->handle, INFINITE);
    switch (result) {
    case WAIT_OBJECT_0:
        if (value_ptr) *value_ptr = thread->arg;
        return 0;
    case WAIT_ABANDONED:
        return EINVAL;
    default:
        return GetLastError();
    }
}

#ifdef ZSTD_CUSTOM_CV
/*******************************************************************************
 * The custom condition variable code is derived from:
 *   https://chromium.googlesource.com/chromium/chromium/+/master/base/synchronization/condition_variable_win.cc
 * The implementation details differ slightly to be more idiomatic C.
 *******************************************************************************
 * Copyright (c) 2011 The Chromium Authors. All rights reserved.
 * Use of this source code is governed by the BSD-style license below.
 *******************************************************************************
 * (c) 2012 The Chromium Authors. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following disclaimer
 * in the documentation and/or other materials provided with the
 * distribution.
 *    * Neither the name of Google Inc. nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * THREAD_event_t is a doubly linked list of Windows events.
 * Each waiting thread gets its own event, and after its use the waiting
 * thread returns it to an event list.
 * Each list has two types of nodes:
 *  1. One list node with handle == 0, which is always the head of the ist.
 *  2. The rest of the nodes are item nodes with handle != 0.
 * Only item nodes can be pushed onto and popped off of the list.
 */

typedef struct THREAD_event_s {
  HANDLE handle;
  struct THREAD_event_s* prev;
  struct THREAD_event_s* next;
} THREAD_event_t;

typedef enum { THREAD_event_container, THREAD_event_event } THREAD_event_type;

/**
 * Validation functions: Used only for sanity check asserts
 */

static DEBUG_ONLY int THREAD_event_validate_links(const THREAD_event_t* events) {
  return events && events->next->prev == events && events->prev->next == events;
}

static DEBUG_ONLY int THREAD_event_validate_list(const THREAD_event_t* events) {
  return events && events->handle == 0 && THREAD_event_validate_links(events);
}

static DEBUG_ONLY int THREAD_event_validate_item(const THREAD_event_t* events) {
  return events && events->handle != 0 && THREAD_event_validate_links(events);
}

static DEBUG_ONLY int THREAD_event_validate_distinct(const THREAD_event_t* event,
                                          const THREAD_event_t* other) {
  return THREAD_event_validate_links(event)
      && THREAD_event_validate_links(other)
      && event != other;
}

/**
 * Does this `THREAD_event_t` only contain item?
 * (It could either be a list or an item)
 */
static int THREAD_event_singleton(const THREAD_event_t* events) {
  assert(THREAD_event_validate_links(events));
  return events == events->next;
}

/**
 * Is this an empty list?
 */
static int THREAD_event_empty(const THREAD_event_t* events) {
  assert(THREAD_event_validate_list(events));
  return THREAD_event_singleton(events);
}

/**
 * Remove this *item* from the list (it must not be the container).
 */
static THREAD_event_t* THREAD_event_extract(THREAD_event_t* events) {
  assert(THREAD_event_validate_item(events));
  if (!THREAD_event_singleton(events)) {
    /* stitch neighbors together */
    events->prev->next = events->next;
    events->next->prev = events->prev;
    /* Make extractee into a singleton */
    events->prev = events;
    events->next = events;
  }
  assert(THREAD_event_singleton(events));
  return events;
}

static THREAD_event_t* THREAD_event_pop_front(THREAD_event_t* events) {
  assert(!THREAD_event_empty(events));
  return THREAD_event_extract(events->next);
}

static THREAD_event_t* THREAD_event_pop_back(THREAD_event_t* events) {
  assert(!THREAD_event_empty(events));
  return THREAD_event_extract(events->prev);
}

static void THREAD_event_push_back(THREAD_event_t* events,
                                   THREAD_event_t* other) {
  assert(THREAD_event_validate_list(events));
  assert(THREAD_event_singleton(other) && THREAD_event_validate_item(other));
  /* Prepare event for insertion */
  other->prev = events->prev;
  other->next = events;
  /* Cut into list */
  events->prev->next = other;
  events->prev = other;
  assert(THREAD_event_validate_distinct(events, other));
}

/**
 * Create an event item or events container based on type.
 */
static THREAD_event_t* THREAD_event_create(THREAD_event_type type) {
  THREAD_event_t* events = malloc(sizeof(THREAD_event_t));
  if (!events) { return NULL; }
  if (type == THREAD_event_event) {
    /* This is not the head of the list, create an event */
    events->handle = CreateEvent(NULL, 0, 0, NULL);
    assert(events->handle);
  } else {
    /* This is the head of the list */
    assert(type == THREAD_event_container);
    events->handle = 0;
  }
  events->prev = events;
  events->next = events;
  return events;
}

static void THREAD_event_destroy(THREAD_event_t* events) {
  if (!events) { return; }
  if (events->handle == 0) {
    /* This is the head of the list destroy all the items */
    while (!THREAD_event_empty(events)) {
      THREAD_event_t* event = THREAD_event_pop_front(events);
      THREAD_event_destroy(event);
    }
  }
  /* It is now either an empty list or an item so only destory
   * the one element.
   */
  assert(THREAD_event_singleton(events));
  if (events->handle != 0) {
    CloseHandle(events->handle);
  }
  free(events);
}

/**
 * Take an event from `from`, put it in `to`, and return it.
 * If `from` is empty create a new event.
 */
static THREAD_event_t* THREAD_event_transfer(THREAD_event_t* from,
                                             THREAD_event_t* to) {
  THREAD_event_t* event;
  if (THREAD_event_empty(from)) {
    event = THREAD_event_create(THREAD_event_event);
    assert(THREAD_event_validate_item(event));
  } else {
    /* FIFO for recycled events to minimize spurious wakeups */
    event = THREAD_event_pop_front(from);
  }
  THREAD_event_push_back(to, event);
  return event;
}

int _pthread_cond_init(pthread_cond_t* cond, const void* attrs) {
  (void)attrs;
  if (!cond) { return 1; }
  pthread_mutex_init(&cond->internalMutex, NULL);
  cond->waitingEvents = THREAD_event_create(THREAD_event_container);
  cond->recycledEvents = THREAD_event_create(THREAD_event_container);
  if (!cond->waitingEvents || !cond->recycledEvents) {
    _pthread_cond_destroy(cond);
    return 1;
  }
  cond->state = RUNNING;
  return 0;
}

int _pthread_cond_destroy(pthread_cond_t* cond) {
  if (!cond) { return 1; }
  /* Disallow any further waiting on this condition variable */
  cond->state = SHUTDOWN;
  pthread_mutex_destroy(&cond->internalMutex);
  THREAD_event_destroy(cond->waitingEvents);
  THREAD_event_destroy(cond->recycledEvents);
  return 0;
}

/**
 * Grab an event (either recycled or new).
 * Push the event onto the waitingEvents list.
 * Release the externalMutex.
 * Wait on the event.
 * Recycle the event.
 * Reacquire the externalMutex.
 */
int _pthread_cond_wait(pthread_cond_t* cond, pthread_mutex_t* externalMutex) {
  THREAD_event_t* waitingEvent = NULL;
  HANDLE handle = 0;
  int bad = 0;

  /* Get an event together with its handle.
   * Errors are delayed until after the critical section.
   */
  pthread_mutex_lock(&cond->internalMutex);
  {
    /* Check if cond is being destroyed or is already destroyed */
    bad |= cond->state != RUNNING;
    if (!bad) {
      /* Grab an event and put it in the waitingEvents list */
      waitingEvent = THREAD_event_transfer(cond->recycledEvents,
                                           cond->waitingEvents);
      bad |= !waitingEvent;
    }
    if (!bad) {
      /* Check the handle to make sure it is valid */
      handle = waitingEvent->handle;
      bad |= !handle;
    }
  }
  pthread_mutex_unlock(&cond->internalMutex);
  if (bad) { return bad; }

  /* Unlock the external mutex and wait for the event to be triggered. */
  pthread_mutex_unlock(externalMutex);
  {
    WaitForSingleObject(handle, INFINITE);
    /* Recycle event immediately to minimize number of events created */
    pthread_mutex_lock(&cond->internalMutex);
    {
      THREAD_event_push_back(cond->recycledEvents, waitingEvent);
    }
    pthread_mutex_unlock(&cond->internalMutex);
  }
  pthread_mutex_lock(externalMutex);
  return 0;
}

/**
 * If there are no waiting threads then do nothing.
 * Otherwise take one threads event off of waitingEvents and signal it.
 */
int _pthread_cond_signal(pthread_cond_t* cond) {
  HANDLE handle;

  pthread_mutex_lock(&cond->internalMutex);
  if (THREAD_event_empty(cond->waitingEvents)) {
    /* No threads are waiting => no op */
    handle = 0;
  } else {
    /* LIFO for waiting events.
     * We take the event for the thread that has most recently
     * waited on the condition variable and signal it.
     *
     * The event is returned to the recycledEvents list by the
     * waiting thread, so we can forget about it.
     */
    handle = THREAD_event_pop_back(cond->waitingEvents)->handle;
    assert(handle);
  }
  pthread_mutex_unlock(&cond->internalMutex);
  /* Signal the event if it exists */
  if (handle) {
    SetEvent(handle);
  }
  return 0;
}

/**
 * Create an empty events list.
 * Swap the waitingEvents list wit the empty events list.
 * For each
 */
int _pthread_cond_broadcast(pthread_cond_t* cond) {
  /* Create an empty events list */
  THREAD_event_t* events = THREAD_event_create(THREAD_event_container);
  if (!events) { return 1; }
  /* Swap the empty events list with the waiting events */
  pthread_mutex_lock(&cond->internalMutex);
  if (!THREAD_event_empty(cond->waitingEvents)) {
    THREAD_event_t* const tmp = cond->waitingEvents;
    cond->waitingEvents = events;
    events = tmp;
  }
  pthread_mutex_unlock(&cond->internalMutex);
  /* Signal every waiting event */
  while (!THREAD_event_empty(events)) {
    /* LIFO for waiting events.
     * Take the most recent event and signal it.
     * Repeat until all waiting events have been signaled.
     *
     * The event is returned to the recycledEvents list by the
     * waiting thread, so we can forget about it.
     */
    HANDLE handle = THREAD_event_pop_back(events)->handle;
    SetEvent(handle);
  }
  /* Destroy the now empty local list we created */
  THREAD_event_destroy(events);
  return 0;
}

#endif

#endif
