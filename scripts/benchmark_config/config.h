#ifndef MESSI_BENCHMARK_CONFIG_H
#define MESSI_BENCHMARK_CONFIG_H

/* The standalone SIMD benchmark does not require a configured MESSI build.
 * ADS_HAVE_AVX2 is supplied explicitly by its compiler command. */
#ifndef ADS_HAVE_AVX2
#define ADS_HAVE_AVX2 0
#endif

#ifndef HAVE_PTHREAD_BARRIER
#define HAVE_PTHREAD_BARRIER 1
#endif

#endif
