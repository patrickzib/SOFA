#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ads/lower_bound_simd.h"
#include "ads/sax/ts.h"

#define DEFAULT_CANDIDATES 2000000ULL
#define DEFAULT_TRIALS 5
#define MAX_DIMS 64
#define MAX_LENGTH 256

typedef struct {
    size_t candidates;
    int threads;
    int trials;
    uint32_t seed;
    int random_access;
} benchmark_options;

static float record_table[MESSI_RECORD_LB_MAX_DIMENSIONS][256];
static sax_type *words;
static ts_type *series;
static ts_type query[MAX_LENGTH];
static uint32_t *order;
static volatile float result_sink;
static uint32_t random_state;

static uint32_t next_random(void) {
    uint32_t value = random_state;
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    random_state = value;
    return value;
}

static double now_seconds(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double) now.tv_sec + (double) now.tv_nsec * 1e-9;
}

static int parse_size(const char *text, size_t *result) {
    char *end = NULL;
    errno = 0;
    unsigned long long value = strtoull(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value == 0 ||
        value > UINT32_MAX || value > SIZE_MAX) return 0;
    *result = (size_t) value;
    return 1;
}

static int parse_positive_int(const char *text, int *result) {
    char *end = NULL;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || value <= 0 || value > INT_MAX)
        return 0;
    *result = (int) value;
    return 1;
}

static void print_usage(const char *program) {
    fprintf(stderr,
        "Usage: %s [OPTIONS]\n"
        "  --count N       Candidate records (default: 2000000)\n"
        "  --threads N     OpenMP workers (default: 1)\n"
        "  --trials N      Timed trials per kernel; median is reported (default: 5)\n"
        "  --seed N        Deterministic 32-bit random seed (default: 1)\n"
        "  --sequential    Scan records in allocation order instead of randomized order\n"
        "  --help          Show this help\n",
        program);
}

static int parse_options(int argc, char **argv, benchmark_options *options) {
    *options = (benchmark_options) {
        .candidates = (size_t) DEFAULT_CANDIDATES,
        .threads = 1,
        .trials = DEFAULT_TRIALS,
        .seed = 1,
        .random_access = 1
    };
    for (int argument = 1; argument < argc; ++argument) {
        if (strcmp(argv[argument], "--help") == 0) {
            print_usage(argv[0]);
            exit(0);
        } else if (strcmp(argv[argument], "--sequential") == 0) {
            options->random_access = 0;
        } else if (strcmp(argv[argument], "--count") == 0 && argument + 1 < argc) {
            if (!parse_size(argv[++argument], &options->candidates)) return 0;
        } else if (strcmp(argv[argument], "--threads") == 0 && argument + 1 < argc) {
            if (!parse_positive_int(argv[++argument], &options->threads)) return 0;
        } else if (strcmp(argv[argument], "--trials") == 0 && argument + 1 < argc) {
            if (!parse_positive_int(argv[++argument], &options->trials)) return 0;
        } else if (strcmp(argv[argument], "--seed") == 0 && argument + 1 < argc) {
            size_t seed = 0;
            if (!parse_size(argv[++argument], &seed)) return 0;
            options->seed = (uint32_t) seed;
        } else {
            return 0;
        }
    }
    return 1;
}

static int allocate_inputs(size_t candidates) {
    if (candidates > SIZE_MAX / MAX_DIMS ||
        candidates * MAX_DIMS > SIZE_MAX / sizeof(*words) ||
        candidates > SIZE_MAX / MAX_LENGTH ||
        candidates * MAX_LENGTH > SIZE_MAX / sizeof(*series) ||
        candidates > SIZE_MAX / sizeof(*order)) return 0;
    return posix_memalign((void **) &words, 64,
                         candidates * MAX_DIMS * sizeof(*words)) == 0 &&
           posix_memalign((void **) &series, 64,
                         candidates * MAX_LENGTH * sizeof(*series)) == 0 &&
           posix_memalign((void **) &order, 64,
                         candidates * sizeof(*order)) == 0;
}

static void initialize_inputs(const benchmark_options *options) {
    random_state = options->seed;
    for (int dimension = 0; dimension < MAX_DIMS; ++dimension) {
        for (int symbol = 0; symbol < 256; ++symbol) {
            record_table[dimension][symbol] =
                (float) (next_random() & 0xffffU) * (1.0f / 65536.0f);
        }
    }
    for (size_t value = 0; value < options->candidates * MAX_DIMS; ++value)
        words[value] = (sax_type) next_random();
    for (size_t value = 0; value < options->candidates * MAX_LENGTH; ++value)
        series[value] = (float) (int32_t) next_random() * (1.0f / 2147483648.0f);
    for (int value = 0; value < MAX_LENGTH; ++value)
        query[value] = (float) (int32_t) next_random() * (1.0f / 2147483648.0f);
    for (size_t candidate = 0; candidate < options->candidates; ++candidate)
        order[candidate] = (uint32_t) candidate;
    if (options->random_access) {
        for (size_t candidate = options->candidates - 1; candidate > 0; --candidate) {
            size_t selected = (size_t) next_random() % (candidate + 1);
            uint32_t temporary = order[candidate];
            order[candidate] = order[selected];
            order[selected] = temporary;
        }
    }
}

static double run_lower_bound(const benchmark_options *options, int dimensions) {
    float sum = 0.0f;
    const double start = now_seconds();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(options->threads) reduction(+:sum)
#endif
    for (size_t candidate = 0; candidate < options->candidates; ++candidate) {
        const uint32_t record = order[candidate];
        sum += messi_record_lb_table_sum(record_table,
            words + (size_t) record * MAX_DIMS, dimensions, FLT_MAX, 1);
    }
    const double elapsed = now_seconds() - start;
    result_sink = sum;
    return elapsed;
}

static double run_exact_distance(const benchmark_options *options, int length) {
    float sum = 0.0f;
    const double start = now_seconds();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(options->threads) reduction(+:sum)
#endif
    for (size_t candidate = 0; candidate < options->candidates; ++candidate) {
        const uint32_t record = order[candidate];
        sum += ts_euclidean_distance_SIMD(query,
            series + (size_t) record * length, length, FLT_MAX);
    }
    const double elapsed = now_seconds() - start;
    result_sink = sum;
    return elapsed;
}

static int compare_double(const void *left, const void *right) {
    const double a = *(const double *) left;
    const double b = *(const double *) right;
    return (a > b) - (a < b);
}

static double median_lower_bound(const benchmark_options *options, int dimensions) {
    double *samples = malloc((size_t) options->trials * sizeof(*samples));
    if (samples == NULL) return -1.0;
    for (int trial = 0; trial < options->trials; ++trial)
        samples[trial] = run_lower_bound(options, dimensions);
    qsort(samples, (size_t) options->trials, sizeof(*samples), compare_double);
    const double median = samples[options->trials / 2];
    free(samples);
    return median;
}

static double median_exact_distance(const benchmark_options *options, int length) {
    double *samples = malloc((size_t) options->trials * sizeof(*samples));
    if (samples == NULL) return -1.0;
    for (int trial = 0; trial < options->trials; ++trial)
        samples[trial] = run_exact_distance(options, length);
    qsort(samples, (size_t) options->trials, sizeof(*samples), compare_double);
    const double median = samples[options->trials / 2];
    free(samples);
    return median;
}

static const char *simd_backend(void) {
#if defined(__AVX512F__)
    return "AVX-512";
#elif ADS_HAVE_AVX2
    return "AVX2";
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
    return "NEON";
#else
    return "scalar";
#endif
}

int main(int argc, char **argv) {
    benchmark_options options;
    if (!parse_options(argc, argv, &options)) {
        print_usage(argv[0]);
        return 2;
    }
#ifndef _OPENMP
    if (options.threads != 1) {
        fputs("This binary was built without OpenMP; --threads must be 1.\n", stderr);
        return 2;
    }
#endif
    if (!allocate_inputs(options.candidates)) {
        fprintf(stderr, "Unable to allocate %.2f GiB of benchmark data.\n",
            ((double) options.candidates *
             (MAX_DIMS * sizeof(*words) + MAX_LENGTH * sizeof(*series) + sizeof(*order))) /
                (1024.0 * 1024.0 * 1024.0));
        free(order);
        free(series);
        free(words);
        return 1;
    }
    initialize_inputs(&options);

    /* Start the OpenMP worker pool and warm the instruction paths without
     * consuming another complete multi-million-record trial. */
    benchmark_options warmup = options;
    if (warmup.candidates > 65536) warmup.candidates = 65536;
    (void) run_lower_bound(&warmup, MAX_DIMS);
    (void) run_exact_distance(&warmup, MAX_LENGTH);

    printf("backend=%s candidates=%zu threads=%d trials=%d access=%s "
           "early_abandon=off data=%.2f_GiB\n",
           simd_backend(), options.candidates, options.threads, options.trials,
           options.random_access ? "random" : "sequential",
           ((double) options.candidates *
            (MAX_DIMS * sizeof(*words) + MAX_LENGTH * sizeof(*series) + sizeof(*order))) /
               (1024.0 * 1024.0 * 1024.0));
    puts("kernel,size,ns_per_candidate,million_candidates_per_second");

    const int dimensions[] = {16, 32, 48, 64};
    double lower_bound_seconds[sizeof(dimensions) / sizeof(dimensions[0])];
    for (size_t width = 0; width < sizeof(dimensions) / sizeof(dimensions[0]); ++width) {
        lower_bound_seconds[width] = median_lower_bound(&options, dimensions[width]);
        if (lower_bound_seconds[width] < 0.0) return 1;
        const double per_candidate = lower_bound_seconds[width] / options.candidates;
        printf("LB,%d,%.3f,%.3f\n", dimensions[width], per_candidate * 1e9,
               1e-6 / per_candidate);
    }

    const int lengths[] = {100, 256};
    double exact_seconds[sizeof(lengths) / sizeof(lengths[0])];
    for (size_t length = 0; length < sizeof(lengths) / sizeof(lengths[0]); ++length) {
        exact_seconds[length] = median_exact_distance(&options, lengths[length]);
        if (exact_seconds[length] < 0.0) return 1;
        const double per_candidate = exact_seconds[length] / options.candidates;
        printf("ED,%d,%.3f,%.3f\n", lengths[length], per_candidate * 1e9,
               1e-6 / per_candidate);
    }

    puts("break_even_lb_dims,ed_length,minimum_percent_pruned");
    for (size_t width = 0; width < sizeof(dimensions) / sizeof(dimensions[0]); ++width) {
        for (size_t length = 0; length < sizeof(lengths) / sizeof(lengths[0]); ++length) {
            printf("%d,%d,%.2f\n", dimensions[width], lengths[length],
                   100.0 * lower_bound_seconds[width] / exact_seconds[length]);
        }
    }

    free(order);
    free(series);
    free(words);
    return result_sink == -1.0f;
}
