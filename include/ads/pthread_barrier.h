#ifndef ADS_PTHREAD_BARRIER_COMPAT_H
#define ADS_PTHREAD_BARRIER_COMPAT_H

#include <pthread.h>
#include <errno.h>

#ifndef HAVE_PTHREAD_BARRIER
#define HAVE_PTHREAD_BARRIER 0
#endif

#if !HAVE_PTHREAD_BARRIER

#ifndef PTHREAD_BARRIER_SERIAL_THREAD
#define PTHREAD_BARRIER_SERIAL_THREAD 1
#endif

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    unsigned count;
    unsigned trip_count;
    unsigned generation;
} pthread_barrier_t;

typedef int pthread_barrierattr_t;

static inline int pthread_barrier_init(pthread_barrier_t *barrier,
                                       const pthread_barrierattr_t *attr,
                                       unsigned count) {
    (void) attr;
    if (barrier == NULL || count == 0)
        return EINVAL;
    barrier->count = 0;
    barrier->trip_count = count;
    barrier->generation = 0;
    pthread_mutex_init(&barrier->mutex, NULL);
    pthread_cond_init(&barrier->cond, NULL);
    return 0;
}

static inline int pthread_barrier_destroy(pthread_barrier_t *barrier) {
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

static inline int pthread_barrier_wait(pthread_barrier_t *barrier) {
    pthread_mutex_lock(&barrier->mutex);
    unsigned generation = barrier->generation;

    barrier->count++;
    if (barrier->count == barrier->trip_count) {
        barrier->generation++;
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return PTHREAD_BARRIER_SERIAL_THREAD;
    }

    while (generation == barrier->generation) {
        pthread_cond_wait(&barrier->cond, &barrier->mutex);
    }
    pthread_mutex_unlock(&barrier->mutex);
    return 0;
}

#endif /* !HAVE_PTHREAD_BARRIER */

#endif /* ADS_PTHREAD_BARRIER_COMPAT_H */
