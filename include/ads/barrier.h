#ifndef ADS_BARRIER_H
#define ADS_BARRIER_H

#include <pthread.h>

typedef struct ads_barrier {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    unsigned int count;
    unsigned int threshold;
    unsigned int generation;
} ads_barrier_t;

static inline int ads_barrier_init(ads_barrier_t *barrier, unsigned int threshold) {
    if (barrier == NULL || threshold == 0) {
        return -1;
    }
    pthread_mutex_init(&barrier->mutex, NULL);
    pthread_cond_init(&barrier->cond, NULL);
    barrier->count = 0;
    barrier->threshold = threshold;
    barrier->generation = 0;
    return 0;
}

static inline int ads_barrier_destroy(ads_barrier_t *barrier) {
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

static inline int ads_barrier_wait(ads_barrier_t *barrier) {
    pthread_mutex_lock(&barrier->mutex);
    unsigned int generation = barrier->generation;

    barrier->count++;
    if (barrier->count == barrier->threshold) {
        barrier->generation++;
        barrier->count = 0;
        pthread_cond_broadcast(&barrier->cond);
        pthread_mutex_unlock(&barrier->mutex);
        return 1;
    }

    while (generation == barrier->generation) {
        pthread_cond_wait(&barrier->cond, &barrier->mutex);
    }
    pthread_mutex_unlock(&barrier->mutex);
    return 0;
}

#endif
