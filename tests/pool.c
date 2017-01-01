#include "pool.h"
#include "threading.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define ASSERT_TRUE(p)                                                         \
  do {                                                                         \
    if (!(p)) {                                                                \
      return 1;                                                                \
    }                                                                          \
  } while (0)
#define ASSERT_FALSE(p) ASSERT_TRUE(!(p))
#define ASSERT_EQ(lhs, rhs) ASSERT_TRUE((lhs) == (rhs))

struct data {
  pthread_mutex_t mutex;
  unsigned data[128];
  size_t i;
};

void fn(void *opaque) {
  struct data *data = (struct data *)opaque;
  if (pthread_mutex_lock(&data->mutex)) {
    exit(2);
  }
  data->data[data->i] = data->i;
  ++data->i;
  if (pthread_mutex_unlock(&data->mutex)) {
    exit(3);
  }
}

int testOrder(size_t numThreads, size_t queueLog) {
  struct data data;
  POOL_ctx *ctx = POOL_create(numThreads, queueLog);
  ASSERT_TRUE(ctx);
  data.i = 0;
  if (pthread_mutex_init(&data.mutex, NULL)) {
    exit(4);
  }
  {
    size_t i;
    for (i = 0; i < 128; ++i) {
      POOL_add(ctx, &fn, &data);
    }
  }
  POOL_free(ctx);
  ASSERT_EQ(128, data.i);
  {
    size_t i;
    for (i = 0; i < data.i; ++i) {
      ASSERT_EQ(i, data.data[i]);
    }
  }
  if (pthread_mutex_destroy(&data.mutex)) {
    exit(5);
  }
  return 0;
}

int main(int argc, const char **argv) {
  size_t numThreads;
  for (numThreads = 1; numThreads <= 4; ++numThreads) {
    size_t queueSize;
    for (queueSize = 1; queueSize <= 2; ++queueSize) {
      if (testOrder(numThreads, queueSize)) {
        printf("FAIL: testOrder\n");
        return 1;
      }
    }
  }
  printf("PASS: testOrder\n");
  (void)argc;
  (void)argv;
  return (POOL_create(0, 1) || POOL_create(1, 0))
      ? printf("FAIL: testInvalid\n"), 1
      : printf("PASS: testInvalid\n"), 0;
}
