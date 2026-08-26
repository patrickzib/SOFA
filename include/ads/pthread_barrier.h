#ifndef ADS_PTHREAD_BARRIER_COMPAT_H
#define ADS_PTHREAD_BARRIER_COMPAT_H

#if !(defined(HAVE_PTHREAD_BARRIER) && HAVE_PTHREAD_BARRIER)

#include "ads/barrier.h"

#ifndef PTHREAD_BARRIER_SERIAL_THREAD
#define PTHREAD_BARRIER_SERIAL_THREAD 1
#endif

#undef pthread_barrier_t
#undef pthread_barrier_init
#undef pthread_barrier_wait
#undef pthread_barrier_destroy

#define pthread_barrier_t ads_barrier_t
#define pthread_barrier_init(barrier, attr, count) ads_barrier_init((barrier), (unsigned int)(count))
#define pthread_barrier_destroy(barrier) ads_barrier_destroy((barrier))
#define pthread_barrier_wait(barrier) ads_barrier_wait((barrier))

#endif /* !(HAVE_PTHREAD_BARRIER) */

#endif /* ADS_PTHREAD_BARRIER_COMPAT_H */
