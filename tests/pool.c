#include "pool.h"
#include "threading.h"
#include "util.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define ASSERT_TRUE(p)                                                         \
  do {                                                                         \
    if (!(p)) {                                                                \
      printf("Assert failed on line: %d\n", __LINE__);                         \
      return 1;                                                                \
    }                                                                          \
  } while (0)
#define ASSERT_FALSE(p) ASSERT_TRUE(!(p))
#define ASSERT_EQ(lhs, rhs) ASSERT_TRUE((lhs) == (rhs))
#define ASSERT_GE(lhs, rhs) ASSERT_TRUE((lhs) >= (rhs))
#define ASSERT_LE(lhs, rhs) ASSERT_TRUE((lhs) <= (rhs))
#define ATOMIC(mutex, macro)                                                   \
  do {                                                                         \
    pthread_mutex_lock((mutex));                                               \
    macro;                                                                     \
    pthread_mutex_unlock((mutex));                                             \
  } while (0)

#define RUN(test, ...)                                                         \
  do {                                                                         \
    int result = test(__VA_ARGS__);                                            \
    if (result) {                                                              \
      printf("FAIL: %s\n", #test);                                             \
      return result ;                                                          \
    } else {                                                                   \
      printf("PASS: %s\n", #test);                                             \
    }                                                                          \
  } while (0)

typedef struct {
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  unsigned data[16];
  size_t i;
  size_t size;
} data_t;

void data_init(data_t *data) {
  pthread_mutex_init(&data->mutex, NULL);
  pthread_cond_init(&data->cond, NULL);
  data->i = 0;
  data->size = sizeof(data->data) / sizeof(*data->data);
}

void data_destroy(data_t *data) {
  pthread_mutex_destroy(&data->mutex);
  pthread_cond_destroy(&data->cond);
}

void update(void *opaque) {
  data_t *data = (data_t *)opaque;
  pthread_mutex_lock(&data->mutex);
  data->data[data->i] = data->i;
  ++data->i;
  pthread_mutex_unlock(&data->mutex);
}

int testOrder(size_t numThreads, size_t queueSize) {
  data_t data;
  POOL_ctx *ctx = POOL_create(numThreads, queueSize);
  ASSERT_TRUE(ctx);
  data_init(&data);
  {
    size_t i;
    for (i = 0; i < data.size; ++i) {
      POOL_add(ctx, &update, &data);
    }
  }
  POOL_free(ctx);
  ASSERT_EQ(data.size, data.i);
  {
    size_t i;
    for (i = 0; i < data.i; ++i) {
      ASSERT_EQ(i, data.data[i]);
    }
  }
  data_destroy(&data);
  return 0;
}

static int testInvalid(void) {
  ASSERT_FALSE(POOL_create(0, 1));
  ASSERT_FALSE(POOL_create(1, 0));
  return 0;
}

static void waiter(void *opaque) {
  data_t* data = (data_t*)opaque;
  pthread_mutex_lock(&data->mutex);
  ++data->data[1];
  while (data->data[0] == 0) {
    pthread_cond_wait(&data->cond, &data->mutex);
  }
  --data->data[0];
  ++data->data[2];
  pthread_mutex_unlock(&data->mutex);
}

static void add_jobs(data_t *data, unsigned jobs, unsigned threads) {
  pthread_mutex_lock(&data->mutex);
  data->data[0] += jobs;
  while (data->data[1] < threads) {
    pthread_mutex_unlock(&data->mutex);
    UTIL_sleepMilli(30);
    pthread_mutex_lock(&data->mutex);
  }
  data->data[1] -= jobs;
  pthread_mutex_unlock(&data->mutex);
}

static int testWaitBroadcast(void) {
  POOL_ctx *ctx;
  data_t data;
  data_init(&data);
  data.data[0] = 0;
  data.data[1] = 0;
  data.data[2] = 0;
  ctx = POOL_create(4, 10);
  ASSERT_TRUE(ctx);

  {
    size_t i;
    for (i = 0; i < 14; ++i) {
      POOL_add(ctx, &waiter, &data);
    }
  }
  add_jobs(&data, 4, 4);
  pthread_cond_broadcast(&data.cond);
  pthread_cond_broadcast(&data.cond);

  add_jobs(&data, 4, 4);
  ATOMIC(&data.mutex, ASSERT_GE(data.data[2], 4));
  pthread_cond_broadcast(&data.cond);

  add_jobs(&data, 2, 4);
  ATOMIC(&data.mutex, ASSERT_GE(data.data[2], 8));
  pthread_cond_broadcast(&data.cond);

  add_jobs(&data, 0, 4);
  pthread_cond_broadcast(&data.cond);
  ATOMIC(&data.mutex, ASSERT_EQ(data.data[2], 10));

  add_jobs(&data, 4, 4);
  pthread_cond_broadcast(&data.cond);

  POOL_free(ctx);
  ATOMIC(&data.mutex, ASSERT_EQ(data.data[2], 14));
  data_destroy(&data);
  return 0;
}

static int testWaitSignal(void) {
  POOL_ctx *ctx;
  data_t data;
  data_init(&data);
  data.data[0] = 0;
  data.data[1] = 0;
  data.data[2] = 0;
  ctx = POOL_create(4, 8);
  ASSERT_TRUE(ctx);

  {
    size_t i;
    for (i = 0; i < 8; ++i) {
      POOL_add(ctx, &waiter, &data);
    }
  }
  add_jobs(&data, 1, 1);
  pthread_cond_signal(&data.cond);
  add_jobs(&data, 0, 4);
  ATOMIC(&data.mutex, ASSERT_EQ(data.data[2], 1));

  add_jobs(&data, 3, 4);
  pthread_cond_signal(&data.cond);
  pthread_cond_signal(&data.cond);
  pthread_cond_signal(&data.cond);
  pthread_cond_signal(&data.cond);
  add_jobs(&data, 0, 4);
  ATOMIC(&data.mutex, ASSERT_EQ(data.data[2], 4));

  add_jobs(&data, 4, 4);
  ATOMIC(&data.mutex, ASSERT_EQ(data.data[2], 4));
  pthread_cond_signal(&data.cond);
  ATOMIC(&data.mutex, ASSERT_LE(data.data[2], 5));
  pthread_cond_signal(&data.cond);
  ATOMIC(&data.mutex, ASSERT_LE(data.data[2], 6));
  pthread_cond_broadcast(&data.cond);

  POOL_free(ctx);
  ATOMIC(&data.mutex, ASSERT_EQ(data.data[2], 8));
  data_destroy(&data);
  return 0;
}

int main(int argc, const char **argv) {
  size_t numThreads;
  for (numThreads = 1; numThreads <= 4; ++numThreads) {
    size_t queueSize;
    for (queueSize = 1; queueSize <= 2; ++queueSize) {
      RUN(testOrder, numThreads, queueSize);
    }
  }
  RUN(testInvalid);
  RUN(testWaitBroadcast);
  RUN(testWaitSignal);
  (void)argc;
  (void)argv;
  return 0;
}
