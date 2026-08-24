//
//  MESSI CLI entrypoint
//
//  Created by Botao Peng, March 2020.
//  
//
#define _GNU_SOURCE

#define PRODUCT "----------------------------------------------\
\nThis is the Adaptive Leaf iSAX index.\n\
----------------------------------------------\n\n"
#ifdef VALUES
#include <values.h>
#endif

#include "config.h"
#include "../../globals.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#ifdef __linux__
#include <sched.h>
#include <glob.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>

#include "ads/sax/sax.h"
#include "ads/sax/ts.h"
#include "ads/isax_visualize_index.h"
#include "ads/isax_file_loaders.h"
#include "ads/isax_visualize_index.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/isax_query_engine.h"
#include "ads/inmemory_query_engine.h"
#include "ads/parallel_inmemory_query_engine.h"
#include "ads/inmemory_index_engine.h"
#include "ads/parallel_query_engine.h"
#include "ads/parallel_index_engine.h"
#include "ads/inmemory_topk_engine.h"
#include "ads/DTWfunction.h"
#include "ads/sfa/dft.h"
#include "ads/sfa/sfa.h"
#include "ads/spartan/spartan.h"
#include "ads/pisa/pisa.h"
#include "ads/trie/trie.h"
#include "include/ads/isax_file_loaders.h"
//#define PROGRESS_CALCULATE_THREAD_NUMBER 4
//#define PROGRESS_FLUSH_THREAD_NUMBER 4
//#define QUERIES_THREAD_NUMBER 4
//#define DISK_BUFFER_SIZE 32
#define FILENAME_LENGTH 256

isax_index *idx;
int query_report_interval = 10;

void INThandler(int);

static int fanout_to_bits(int fanout) {
    int bits = 0;
    while (fanout > 1 && (fanout & 1) == 0) {
        fanout >>= 1;
        ++bits;
    }
    return fanout == 1 ? bits : -1;
}

static FILE *open_logfile_or_tmp(const char *path) {
    FILE *file = fopen(path, "w");
    if (file != NULL) {
        return file;
    }

    fprintf(stderr, "warning: cannot open log file %s; using a temporary log.\n", path);
    file = tmpfile();
    if (file == NULL) {
        perror("tmpfile");
    }
    return file;
}

static int ensure_directory(const char *path) {
    char directory[FILENAME_LENGTH];
    size_t length = strlen(path);

    if (length == 0 || length >= sizeof(directory)) return -1;
    memcpy(directory, path, length + 1);

    for (char *cursor = directory + 1; *cursor != '\0'; ++cursor) {
        if (*cursor != '/') continue;
        *cursor = '\0';
        if (mkdir(directory, 0777) != 0 && errno != EEXIST) return -1;
        *cursor = '/';
    }

    if (mkdir(directory, 0777) != 0 && errno != EEXIST) return -1;
    return 0;
}

static double monotonic_seconds(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double) now.tv_sec + (double) now.tv_nsec / 1000000000.0;
}

static void format_compact_count(double value, char *buffer, size_t buffer_size) {
    const char *suffix = "";
    double scaled = value;
    if (fabs(value) >= 1000000000.0) { scaled = value / 1000000000.0; suffix = "B"; }
    else if (fabs(value) >= 1000000.0) { scaled = value / 1000000.0; suffix = "M"; }
    else if (fabs(value) >= 1000.0) { scaled = value / 1000.0; suffix = "K"; }
    if (*suffix == '\0') snprintf(buffer, buffer_size, "%.0f", scaled);
    else snprintf(buffer, buffer_size, "%.2f %s", scaled, suffix);
}

/* Keep every query engine, including the trie engines, from discovering a
 * short query file only after an expensive index build. */
static int clamp_query_count_to_file(const char *path, int timeseries_size,
                                     int *query_count) {
    FILE *file;
    long file_size;
    unsigned long long available;
    size_t record_size;

    if (path == NULL || query_count == NULL || *query_count <= 0) return 1;
    if (timeseries_size <= 0) return 0;
    file = fopen(path, "rb");
    if (file == NULL) {
        fprintf(stderr, "error: cannot open query file: %s\n", path);
        return 0;
    }
    if (fseek(file, 0L, SEEK_END) != 0 || (file_size = ftell(file)) < 0) {
        fprintf(stderr, "error: cannot determine size of query file: %s\n", path);
        fclose(file);
        return 0;
    }
    fclose(file);
    record_size = (size_t) timeseries_size * sizeof(ts_type);
    available = (unsigned long long) file_size / record_size;
    if (available == 0) {
        fprintf(stderr, "error: query file contains no complete records: %s\n", path);
        return 0;
    }
    if ((unsigned long long) *query_count > available) {
        fprintf(stderr,
                "warning: query file %s contains only %llu records; reducing --queries-size from %d to %llu.\n",
                path, available, *query_count, available);
        *query_count = available > INT_MAX ? INT_MAX : (int) available;
    }
    return 1;
}

#ifndef __linux__
typedef int cpu_set_t;
static inline void CPU_ZERO(cpu_set_t *set) { (void) set; }
static inline void CPU_SET(int cpu, cpu_set_t *set) { (void) cpu; (void) set; }
static inline int pthread_setaffinity_np(pthread_t thread, size_t size, const cpu_set_t *set) {
    (void) thread;
    (void) size;
    (void) set;
    return 0;
}
static inline int pthread_getaffinity_np(pthread_t thread, size_t size, cpu_set_t *set) {
    (void) thread;
    (void) size;
    (void) set;
    return 0;
}
#endif

#ifdef __linux__
static int parse_cpu_list(const char *text, cpu_set_t *set, const cpu_set_t *allowed) {
    const char *cursor = text;

    while (*cursor != '\0') {
        char *end;
        long first = strtol(cursor, &end, 10);
        long last = first;
        if (end == cursor || first < 0 || first >= CPU_SETSIZE) {
            return -1;
        }
        if (*end == '-') {
            cursor = end + 1;
            last = strtol(cursor, &end, 10);
            if (end == cursor || last < first || last >= CPU_SETSIZE) {
                return -1;
            }
        }
        for (long cpu = first; cpu <= last; ++cpu) {
            if (CPU_ISSET((int) cpu, allowed)) {
                CPU_SET((int) cpu, set);
            }
        }
        if (*end == '\n' || *end == '\0') {
            return 0;
        }
        if (*end != ',') {
            return -1;
        }
        cursor = end + 1;
    }
    return -1;
}

/* Select CPUs belonging to the first requested NUMA nodes that intersect the
 * process cpuset.  A return value of zero leaves the caller's allowed cpuset
 * unchanged, which is the correct fallback on systems without NUMA sysfs. */
static int select_numa_cpus(cpu_set_t *selected, const cpu_set_t *allowed, int requested_nodes) {
    glob_t paths;
    int selected_nodes = 0;
    int result = 0;

    if (glob("/sys/devices/system/node/node[0-9]*", 0, NULL, &paths) != 0) {
        return 0;
    }
    CPU_ZERO(selected);
    for (size_t i = 0; i < paths.gl_pathc &&
                       (requested_nodes < 0 || selected_nodes < requested_nodes); ++i) {
        char filename[PATH_MAX];
        char cpus[4096];
        FILE *file;
        cpu_set_t node_cpus;
        if (snprintf(filename, sizeof(filename), "%s/cpulist", paths.gl_pathv[i]) >=
            (int) sizeof(filename)) {
            continue;
        }
        file = fopen(filename, "r");
        if (file == NULL || fgets(cpus, sizeof(cpus), file) == NULL) {
            if (file != NULL) fclose(file);
            continue;
        }
        fclose(file);
        CPU_ZERO(&node_cpus);
        if (parse_cpu_list(cpus, &node_cpus, allowed) == 0 && CPU_COUNT(&node_cpus) > 0) {
            for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
                if (CPU_ISSET(cpu, &node_cpus)) {
                    CPU_SET(cpu, selected);
                }
            }
            selected_nodes++;
            result = selected_nodes;
        }
    }
    globfree(&paths);
    return result;
}
#endif

int main(int argc, char **argv) {
    typedef enum {
        MESSI_ROOT_SPLIT_DEFAULT,
        MESSI_ROOT_SPLIT_UNIFORM,
        MESSI_ROOT_SPLIT_VARIANCE
    } messi_root_split_mode;
    signal(SIGINT, INThandler);

#ifndef BENCHMARK
    printf(PRODUCT);

#if VERBOSE_LEVEL == 0
    printf("Executing in silent mode. Please wait.\n");
#endif
#endif
    static char index_directory[FILENAME_LENGTH];
    strcpy(index_directory, getenv("HOME"));
    strcat(index_directory, "/index/");

    static char data_directory[FILENAME_LENGTH];
    strcpy(data_directory, getenv("HOME"));
    strcat(data_directory, "/data/data");

    static char query_directory[FILENAME_LENGTH];
    strcpy(query_directory, getenv("HOME"));
    strcat(query_directory, "/data/query");

    static char *dataset = data_directory;
    static char *queries = query_directory;
    static char *index_path = index_directory;
    static char *labelset = index_directory;
    //static int dataset_size = 20100;//1000000;
    //static int dataset_size =   1000000;//simple2
    static long int dataset_size = 6000000;//testbench
    static int queries_size = 10;
    static int time_series_size = 256;
    static int n_segments = 16;
    static int sax_cardinality = 8;
    static int leaf_size = 2000;
    static int min_leaf_size = 10;
    static int initial_lbl_size = 2000;
    static int flush_limit = 200000;
    static int initial_fbl_size = 100;
    static char use_index = 0;
    static int complete_type = 0;
    static int total_loaded_leaves = 1;
    static int tight_bound = 0;
    static int aggressive_check = 0;
    static float minimum_distance = FLT_MAX;
    static int serial_scan = 0;
    static char knnlabel = 0;
    static int min_checked_leaves = -1;
    static int requested_threads = 0;
    static int requested_numa_nodes = -1;
    static char inmemory_flag = 0;
    /* SIMD is enabled automatically when AVX2 support was compiled in.
     * --no-simd remains available for scalar comparison runs. */
    static char SIMD_flag = ADS_HAVE_AVX2 ? 1 : 0;
    static char is_norm = 0;
    static int histogram_type = 1;
    static int sample_type = 1;
    static unsigned int sampling_seed = 1;
    static int n_coefficients = 0;
    static int filetype_int = 0;
    static int apply_znorm = 0;
    static int dynamic_index = 1;
    static messi_root_split_mode root_split_mode = MESSI_ROOT_SPLIT_DEFAULT;
    static int node_split_criterion = 1;
    static messi_index_type index_type = MESSI_INDEX_ISAX;
    /* The trie mirrors iSAX's per-query worker scheduling by default. */
    static int trie_query_batch = 0;
    static int profile_query_phases_requested = 0;
    static int queue_number_specified = 0;
    static int trie_mbr_dimensions = 0;
    static int trie_split_dimensions = 0;
    static int trie_record_mbr_suffix_bound = 0;
    static int trie_fanout = 8;
    static int trie_fanout_specified = 0;
    static int trie_dynamic_alphabet = 0;
    static int trie_min_fanout = 2;
    static int trie_max_fanout = 16;
    static int trie_alphabet_budget_bits = 3;
    static int query_report_interval_requested = 10;

    int calculate_thread = 8;
    int function_type = 0;
    N_PQUEUE = 1;
    maxreadthread = 5;
    read_block_length = 20000;
    int k_size = 0;
    long int labelsize = 1;
    int topk = 0;
    int dtwwindowsize = 0;
    int sample_size = 1000;

    time_t time_now;

    while (1) {
        static struct option long_options[] = {
                {"use-index",           no_argument,       0,    'a'},
                {"initial-lbl-size",    required_argument, 0,    'b'},
                {"complete-type",       required_argument, 0,    'c'},
                {"dataset",             required_argument, 0,    'd'},
                {"total-loaded-leaves", required_argument, 0,    'e'},
                {"flush-limit",         required_argument, 0,    'f'},
                {"aggressive-check",    no_argument,       0,    'g'},
                {"help",                no_argument,       0,    'h'},
                {"initial-fbl-size",    required_argument, 0,    'i'},
                {"serial",              no_argument,       0,    'j'},
                {"queries-size",        required_argument, 0,    'k'},
                {"leaf-size",           required_argument, 0,    'l'},
                {"min-leaf-size",       required_argument, 0,    'm'},
                {"tight-bound",         no_argument,       0,    'n'},
                {"read-thread",         required_argument, 0,    'o'},
                {"index-path",          required_argument, 0,    'p'},
                {"queries",             required_argument, 0,    'q'},
                {"read-block",          required_argument, 0,    'r'},
                {"minimum-distance",    required_argument, 0,    's'},
                {"timeseries-size",     required_argument, 0,    't'},
                {"min-checked-leaves",  required_argument, 0,    'u'},
                {"in-memory",           no_argument,       0,    'v'},
                {"threads",             required_argument, 0,    'I'},
                {"numa",                required_argument, 0,    'J'},
                {"sax-cardinality",     required_argument, 0,    'x'},
                {"n-segments",          required_argument, 0,    'B'},
                {"function-type",       required_argument, 0,    'y'},
                {"dataset-size",        required_argument, 0,    'z'},
                {"k-size",              required_argument, 0,    '0'},
                {"knn-label-set",       required_argument, 0,    '1'},
                {"knn-label-size",      required_argument, 0,    '2'},
                {"knn",                 no_argument,       0,    '3'},
                {"topk",                no_argument,       0,    '4'},
                {"dtwwindowsize",       required_argument, 0,    '5'},
                {"queue-number",        required_argument, 0,    '6'},
                {"sample-size",         required_argument, 0,    '8'},
                {"is-norm",             no_argument,       0,    '9'},
                {"histogram-type",      required_argument, 0,    'A'},
                {"sample-type",         required_argument, 0,    'C'},
                {"sfa-n-coefficients",  required_argument, 0,    'D'},
                {"filetype-int",        no_argument,       0,    'E'},
                {"apply-z-norm",        no_argument,       0,    'F'},
                {"node-split-criterion",required_argument, 0,   'G'},
                {"dynamic-root-split-uniform", required_argument, 0, 'K'},
                {"dynamic-root-split-variance", no_argument,     0, 'L'},
                {"index-type", required_argument, 0, 1000},
                {"trie-query-parallel", no_argument, 0, 1001},
                {"profile-query-phases", no_argument, 0, 1002},
                {"trie-query-batch", no_argument, 0, 1003},
                {"trie-mbr-dimensions", required_argument, 0, 1004},
                {"trie-split-dimensions", required_argument, 0, 1015},
                {"trie-record-mbr-suffix-bound", no_argument, 0, 1013},
                {"trie-fanout", required_argument, 0, 1005},
                {"trie-dynamic-alphabet", no_argument, 0, 1009},
                {"trie-min-fanout", required_argument, 0, 1010},
                {"trie-max-fanout", required_argument, 0, 1011},
                {"trie-alphabet-budget-bits", required_argument, 0, 1012},
                {"query-report-interval", required_argument, 0, 1006},
                {"sampling-seed", required_argument, 0, 1007},
                {"no-simd",             no_argument,       0, 1008},
                {NULL,                  0,                 NULL, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "",
                            long_options, &option_index);
        if (c == -1)
            break;
        switch (c) {
            case 1000:
                if (strcmp(optarg, "isax") == 0) index_type = MESSI_INDEX_ISAX;
                else if (strcmp(optarg, "trie") == 0) index_type = MESSI_INDEX_TRIE;
                else { fprintf(stderr, "error: index-type must be isax or trie.\n"); return EXIT_FAILURE; }
                break;
            case 1001:
                trie_query_batch = 0;
                break;
            case 1002:
                profile_query_phases_requested = 1;
                break;
            case 1003:
                trie_query_batch = 1;
                break;
            case 1004:
                trie_mbr_dimensions = atoi(optarg);
                break;
            case 1015:
                trie_split_dimensions = atoi(optarg);
                break;
            case 1005:
                trie_fanout = atoi(optarg);
                trie_fanout_specified = 1;
                break;
            case 1009:
                trie_dynamic_alphabet = 1;
                break;
            case 1010:
                trie_min_fanout = atoi(optarg);
                break;
            case 1011:
                trie_max_fanout = atoi(optarg);
                break;
            case 1012:
                trie_alphabet_budget_bits = atoi(optarg);
                break;
            case 1006:
                query_report_interval_requested = atoi(optarg);
                break;
            case 1008:
                SIMD_flag = 0;
                break;
            case 1013:
                trie_record_mbr_suffix_bound = 1;
                break;
            case 'j':
                serial_scan = 1;
                break;
            case 'g':
                aggressive_check = 1;
                break;

            case 's':
                minimum_distance = atof(optarg);
                break;

            case 'n':
                tight_bound = 1;
                break;

            case 'e':
                total_loaded_leaves = atoi(optarg);
                break;

            case 'c':
                complete_type = atoi(optarg);
                break;

            case 'q':
                queries = optarg;
                break;

            case 'k':
                queries_size = atoi(optarg);
                break;

            case 'd':
                dataset = optarg;
                break;

            case 'p':
                index_path = optarg;
                break;

            case 'z':
                dataset_size = atoi(optarg);
                break;

            case 't':
                time_series_size = atoi(optarg);
                break;

            case 'x':
                sax_cardinality = atoi(optarg);
                break;

            case 'B':
                n_segments = atoi(optarg);
                break;

            case 'l':
                leaf_size = atoi(optarg);
                break;

            case 'm':
                min_leaf_size = atoi(optarg);
                break;

            case 'b':
                initial_lbl_size = atoi(optarg);
                break;

            case 'f':
                flush_limit = atoi(optarg);
                break;

            case 'u':
                min_checked_leaves = atoi(optarg);
                break;
            case 'I':
                if (strcmp(optarg, "auto") == 0) {
                    requested_threads = 0;
                } else {
                    requested_threads = atoi(optarg);
                    if (requested_threads <= 0) {
                        fprintf(stderr, "error: threads must be a positive integer or 'auto'.\n");
                        return EXIT_FAILURE;
                    }
                }
                break;
            case 'J':
                if (strcmp(optarg, "auto") == 0) {
                    requested_numa_nodes = -1;
                } else if (strcmp(optarg, "none") == 0) {
                    requested_numa_nodes = 0;
                } else {
                    requested_numa_nodes = atoi(optarg);
                    if (requested_numa_nodes <= 0) {
                        fprintf(stderr, "error: numa must be 'auto', 'none', or a positive node count.\n");
                        return EXIT_FAILURE;
                    }
                }
                break;

            case 'y':
                function_type = atoi(optarg);
                break;
            case 'i':
                initial_fbl_size = atoi(optarg);
                break;
            case 'o':
                maxreadthread = atoi(optarg);
                break;
            case 'r':
                read_block_length = atoi(optarg);
                break;
            case '0':
                k_size = atoi(optarg);
                break;
            case '1':
                labelset = optarg;
                break;
            case '2':
                labelsize = atoi(optarg);
                break;
            case '3':
                knnlabel = 1;
                break;
            case '4':
                topk = 1;
                break;
            case '5':
                dtwwindowsize = atoi(optarg);
                break;
            case '6':
                N_PQUEUE = atoi(optarg);
                if (N_PQUEUE <= 0) {
                    fprintf(stderr, "Error: --queue-number must be a positive integer.\n");
                    return EXIT_FAILURE;
                }
                queue_number_specified = 1;
                break;
            case '8':
                sample_size = atoi(optarg);
                break;
            case 'A':
                histogram_type = atoi(optarg);
                break;
            case 'C':
                sample_type = atoi(optarg);
                break;
            case 1007: {
                char *end = NULL;
                unsigned long value = strtoul(optarg, &end, 10);
                if (optarg[0] == '\0' || end == optarg || *end != '\0' ||
                    value > UINT_MAX) {
                    fprintf(stderr, "Error: --sampling-seed must be an unsigned integer.\n");
                    return EXIT_FAILURE;
                }
                sampling_seed = (unsigned int) value;
                break;
            }
            case 'D':
                n_coefficients = atoi(optarg);
                break;
            case 'E':
                filetype_int = 1;
                break;
            case 'F':
                apply_znorm = 1;
                break;
            case 'G':
                node_split_criterion = atoi(optarg);
                if (node_split_criterion < 1 || node_split_criterion > 4) {
                    fprintf(stderr, "error: node-split-criterion must be 1-4 (received %d).\n",
                            node_split_criterion);
                    return EXIT_FAILURE;
                }
                break;
            case 'K':
                if (root_split_mode == MESSI_ROOT_SPLIT_VARIANCE) {
                    fprintf(stderr, "error: --dynamic-root-split-variance cannot be combined with a uniform root split.\n");
                    return EXIT_FAILURE;
                }
                dynamic_index = atoi(optarg);
                root_split_mode = MESSI_ROOT_SPLIT_UNIFORM;
                break;
            case 'L':
                if (root_split_mode == MESSI_ROOT_SPLIT_UNIFORM) {
                    fprintf(stderr, "error: --dynamic-root-split-variance cannot be combined with a uniform root split.\n");
                    return EXIT_FAILURE;
                }
                root_split_mode = MESSI_ROOT_SPLIT_VARIANCE;
                break;

            case 'h':
#ifdef BENCHMARK
                printf(PRODUCT);
#endif
                printf("Usage:\n\
                \t--dataset XX \t\t\tThe path to the dataset file\n\
                \t--queries XX \t\t\tThe path to the queries file\n\
                \t--dataset-size XX \t\tThe number of time series to load\n\
                \t--queries-size XX \t\tThe number of queries to do\n\
                \t--minimum-distance XX\t\tThe minimum distance we search (MAX if not set)\n\
                \t--use-index  \t\t\tSpecifies that an input index will be used\n\
                \t--index-path XX \t\tThe path of the output folder\n\
                \t--timeseries-size XX\t\tThe size of each time series\n\
                \t--sax-cardinality XX\t\tThe maximum sax cardinality in number of bits (power of two).\n\
                \t--n-segments XX\t\tSymbolic dimensions (trie: 16--64; default: 16).\n\
                \t--leaf-size XX\t\t\tThe maximum size of each leaf\n\
                \t--min-leaf-size XX\t\tThe minimum size of each leaf\n\
                \t--initial-lbl-size XX\t\tThe initial lbl buffer size for each buffer.\n\
                \t--flush-limit XX\t\tThe limit of time series in memory at the same time\n\
                \t--initial-fbl-size XX\t\tThe initial fbl buffer size for each buffer.\n\
                \t--node-split-criterion XX\tSelect split decision (1=informed default, 2=simple, 3=maxvar, 4=maxbin)\n\
                \t--dynamic-root-split-uniform XX\tUniform root split in bits per symbolic dimension\n\
                \t--dynamic-root-split-variance\tUse all root-key bits, assigned by transform variance\n\
                \t--index-type isax|trie\tIndex layout (default: isax)\n\
                \t--trie-query-parallel\tParallelize each trie query across subtrees (default)\n\
                \t--trie-query-batch\tBatch independent trie queries instead\n\
                \t--trie-mbr-dimensions XX\tTrie MBR dimensions (16--128; default: maximum available)\n\
                \t--trie-split-dimensions XX\tTrie split-candidate dimensions (default: min(32, MBR dimensions))\n\
                \t--trie-record-mbr-suffix-bound\tAdd non-record-dimension leaf-MBR contributions to trie record bounds\n\
                \t--trie-fanout 2|4|8\tTrie symbolic split fanout (default: 8)\n\
                \t--trie-dynamic-alphabet\tUse one global variance-weighted alphabet allocation\n\
                \t--trie-min-fanout 2|4|...\tMinimum dynamic trie fanout (default: 2)\n\
                \t--trie-max-fanout 2|4|...\tMaximum dynamic trie fanout (default: 16)\n\
                \t--trie-alphabet-budget-bits N\tAverage dynamic alphabet budget in bits (default: 3)\n\
                \t--query-report-interval N\tPrint first, every Nth completed, and final query row (0=none; default: 10)\n\
                \t--profile-query-phases\tRecord direct accumulated worker time for traversal, lower bounds, and exact distances\n\
                \t--complete-type XX\t\t0 for no complete, 1 for serial, 2 for leaf\n\
                \t--total-loaded-leaves XX\tNumber of leaves to load at each fetch\n\
                \t--min-checked-leaves XX\t\tNumber of leaves to check at minimum\n\
                \t--tight-bound XX\t\tSet for tight bounds.\n\
                \t--aggressive-check XX\t\tSet for aggressive check\n\
                \t--serial\t\t\tSet for serial scan\n\
                \t--knn\t\t\t\tSet for knn search\n\
                \t--k-size\t\t\tSet k size of knn\n\
                \t--topk\t\t\t\tSet for top k search\n\
                \t--dtwwindowsize\t\t\tSet dtw window size\n\
                \t--knn-label-set\t\t\tSet label set\n\
                \t--queue-number\t\t\tset the number of priority queues\n\
                \t--in-memory\t\t\tSet for in-memory search\n\
                \t--function-type\t\t\tSet for query answering type\n\
                \t\t\tin memory  only index creation: 0\n\
                \t\t\tParIS-TS: 1\n\
                \t\t\tParIS: 2\n\
                \t\t\t\\MESSI-mq: 3\n\
                \t\t\tMESSI+SFA: 4\n\
                \t\t\tMESSI+SPARTAN: 5\n\
                \t\t\tMESSI+PISA: 6\n\
                \t--no-simd\t\t\tDisable SIMD even when AVX2 is available\n\
                \t--sample-size\t\t\tSet sample size for MCB\n\
                \t--sample-type\t\t\tSet for sampling strategy\n\
                \t\t\tfirst-n-values sampling: 1\n\
                \t\t\tuniform sampling: 2\n\
                \t\t\trandom sampling: 3\n\
                \t--sampling-seed\t\t\tSeed for random binning sampling (default: 1)\n\
                \t--filetype-int\t\t\tSet if the input time series file is stored in int-type\n\
                \t--apply-z-norm\t\t\tApply z-normalization to the data\n\
                \t--is-norm\t\t\tSet for search with normalized input time series\n\
                \t--sfa-n-coefficients\t\t\tSet number of coeff to choose highest-variance coeff (doubled for real & imag parts - must be between n-segments and timeseries-size)\n\
                \t--histogram-type\t\t\tSet for binning strategy\n\
                \t\t\tequi-depth splitting (default): 1\n\
                \t\t\tequi-width splitting: 2\n\
                \t--threads N|auto\t\tWorker threads (default: all CPUs available to this process)\n\
                \t--numa auto|none|N\tCPU affinity: all available NUMA nodes, disabled, or first N nodes\n\
                \t--help\n\n\
                ");
                return 0;
                break;
            case 'a':
                use_index = 1;
                break;
            case 'v':
                inmemory_flag = 1;
                break;
            case '9':
                is_norm = 1;
                break;
            default:
                exit(-1);
                break;
        }
    }
    INIT_STATS();
    profile_query_phases = profile_query_phases_requested;
    if (query_report_interval_requested < 0) {
        fprintf(stderr, "error: query report interval must be zero or positive.\n");
        return EXIT_FAILURE;
    }
    query_report_interval = query_report_interval_requested;

    if (index_type == MESSI_INDEX_TRIE) {
        if (function_type < 3 || function_type > 6) {
            fprintf(stderr, "error: --index-type trie supports function types 3 (SAX), 4 (SFA), 5 (SPARTAN), and 6 (PISA).\n");
            return EXIT_FAILURE;
        }
        if (root_split_mode != MESSI_ROOT_SPLIT_DEFAULT) {
            fprintf(stderr, "error: dynamic root split options are iSAX-only.\n");
            return EXIT_FAILURE;
        }
        if (n_segments < 16 || n_segments > 64) {
            fprintf(stderr, "error: trie n-segments must be between 16 and 64.\n");
            return EXIT_FAILURE;
        }
        if (trie_mbr_dimensions != 0 &&
            (trie_mbr_dimensions < n_segments || trie_mbr_dimensions < 16 ||
             trie_mbr_dimensions > 128 || trie_mbr_dimensions > time_series_size)) {
            fprintf(stderr,
                    "error: trie MBR dimensions must be between n-segments (%d) and min(128, timeseries-size) (%d).\n",
                    n_segments, time_series_size < 128 ? time_series_size : 128);
            return EXIT_FAILURE;
        }
        if (trie_dynamic_alphabet && trie_fanout_specified) {
            fprintf(stderr, "error: --trie-dynamic-alphabet cannot be combined with --trie-fanout.\n");
            return EXIT_FAILURE;
        }
        if (!trie_dynamic_alphabet && trie_fanout != 2 && trie_fanout != 4 && trie_fanout != 8) {
            fprintf(stderr, "error: trie fanout must be 2, 4, or 8.\n");
            return EXIT_FAILURE;
        }
        if (trie_dynamic_alphabet) {
            const int min_bits = fanout_to_bits(trie_min_fanout);
            const int max_bits = fanout_to_bits(trie_max_fanout);
            if (min_bits < 1 || max_bits < min_bits || max_bits > 8) {
                fprintf(stderr, "error: dynamic trie fanouts must be powers of two between 2 and 256, with min <= max.\n");
                return EXIT_FAILURE;
            }
            if (trie_alphabet_budget_bits < min_bits || trie_alphabet_budget_bits > max_bits) {
                fprintf(stderr, "error: dynamic trie average budget must be between min and max fanout bits.\n");
                return EXIT_FAILURE;
            }
        }
    } else if (trie_fanout != 8 || trie_dynamic_alphabet) {
        fprintf(stderr, "error: --trie-fanout requires --index-type trie.\n");
        return EXIT_FAILURE;
    }
    if (trie_record_mbr_suffix_bound && index_type != MESSI_INDEX_TRIE) {
        fprintf(stderr, "error: --trie-record-mbr-suffix-bound requires --index-type trie.\n");
        return EXIT_FAILURE;
    }
    if (dynamic_index < 1 || dynamic_index > sax_cardinality) {
        fprintf(stderr, "error: uniform root split must be between 1 and sax-cardinality (%d).\n",
                sax_cardinality);
        return EXIT_FAILURE;
    }
    if (root_split_mode == MESSI_ROOT_SPLIT_VARIANCE &&
        function_type != 4 && function_type != 5 && function_type != 6) {
        fprintf(stderr, "error: --dynamic-root-split-variance is supported by SFA (4), SPARTAN (5), and PISA (6).\n");
        return EXIT_FAILURE;
    }
    if (!clamp_query_count_to_file(queries, time_series_size, &queries_size)) {
        return EXIT_FAILURE;
    }

#ifdef __linux__
    cpu_set_t allowed_cpus, mask;
    CPU_ZERO(&allowed_cpus);
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(allowed_cpus), &allowed_cpus) != 0 ||
        CPU_COUNT(&allowed_cpus) == 0) {
        perror("sched_getaffinity");
        return EXIT_FAILURE;
    }
    int available_cpus = CPU_COUNT(&allowed_cpus);
    int detected_numa_nodes = 0;
    if (requested_numa_nodes != 0) {
        detected_numa_nodes = select_numa_cpus(&mask, &allowed_cpus, requested_numa_nodes);
        if (detected_numa_nodes == 0 || CPU_COUNT(&mask) == 0) {
            mask = allowed_cpus;
        }
    } else {
        mask = allowed_cpus;
    }
    int usable_cpus = CPU_COUNT(&mask);
    int thread_count = requested_threads > 0 ? requested_threads : usable_cpus;
    if (thread_count > usable_cpus) {
        fprintf(stderr, "warning: requested %d threads but only %d selected CPUs are available; capping threads.\n",
                thread_count, usable_cpus);
        thread_count = usable_cpus;
    }
    if (requested_numa_nodes != 0 && pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask) != 0) {
        perror("pthread_setaffinity_np");
        return EXIT_FAILURE;
    }
    fprintf(stderr, ">>> using %d worker threads on %d available CPUs%s%s\n", thread_count,
            usable_cpus, requested_numa_nodes == 0 ? " (affinity disabled)" : "",
            detected_numa_nodes > 0 ? " across detected NUMA nodes" : "");
#else
    long online_cpus = sysconf(_SC_NPROCESSORS_ONLN);
    int available_cpus = online_cpus > 0 ? (int) online_cpus : 1;
    int thread_count = requested_threads > 0 ? requested_threads : available_cpus;
    if (thread_count > available_cpus) {
        fprintf(stderr, "warning: requested %d threads but only %d CPUs are available; capping threads.\n",
                thread_count, available_cpus);
        thread_count = available_cpus;
    }
    fprintf(stderr, ">>> using %d worker threads; NUMA affinity is unavailable on this platform.\n", thread_count);
#endif
    calculate_thread = thread_count;
    maxquerythread = thread_count;
    /* Keep the default work distribution proportional to the active query
     * workers.  An explicit --queue-number always takes precedence. */
    if (!queue_number_specified) N_PQUEUE = thread_count;
    if (use_index) {
        isax_index *idx = index_read(index_path);
        idx->settings->tight_bound = tight_bound;
        idx->settings->aggressive_check = aggressive_check;
        idx->settings->total_loaded_leaves = total_loaded_leaves;
        idx->settings->min_leaf_size = min_leaf_size;
        print_settings(idx->settings, maxquerythread, trie_query_batch);
        // fprintf(stderr,"total_records: %ld\n", idx->total_records);
        // fprintf(stderr,"loaded_records: %ld\n", idx->loaded_records);
        // create_wedges(idx, NULL);

        if (serial_scan) {
            cache_sax_file(idx);
            if (knnlabel) {
                if (function_type == 0) {
                    isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                               min_checked_leaves, k_size, labelsize,
                                               &exact_topk_serial); //ADS+ knn
                } else if (function_type == 1) {
                    isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                               min_checked_leaves, k_size, labelsize,
                                               &exact_topk_serial_ParIS); //ParIS knn
                }

            } else if (topk) {
                if (function_type == 0) {
                    isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                k_size, &exact_topk_serial);//ADS+ topk
                } else if (function_type == 1) {
                    isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                k_size, &exact_topk_serial_ParIS);//ParIS topk
                }
            } else {
                if (function_type == 0) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search_serial);//ADS+ similarity search
                } else if (function_type == 1) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search_serial_ParIS);//ParIS similarity search
                } else if (function_type == 2) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search_serial_ParIS_nb);//ParIS-nb similarity search
                } else if (function_type == 3) {
                    //isax_query_binary_file_para(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial);
                }
            }

        } else {
            if (knnlabel) {
                if (function_type == 0) {
                    isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                               min_checked_leaves, k_size, labelsize, &exact_topk);
                }

            } else if (topk) {
                if (function_type == 0) {
                    isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                k_size, &exact_topk);
                }
            } else {
                if (function_type == 0) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search);
                } else {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search_m);
                }
            }
        }
        PRINT_STATS(0.00f)

        flush_all_leaf_buffers(idx, TMP_ONLY_CLEAN);
        index_write(idx);

        //clear_wedges(idx, NULL);
        //printf("this is the total search time is %f\n", total_time);
        isax_index_destroy(idx, NULL);
    } else {
        char rm_command[256];
        int index_segments = n_segments;
        int trie_bound_dimensions = 0;

        /* Trie words have extra dimensions available solely for partitioning.
         * The public n-segments value remains the lower-bound prefix. */
        if (index_type == MESSI_INDEX_TRIE) {
            trie_bound_dimensions = index_segments;
            index_segments = trie_mbr_dimensions > 0
                                 ? trie_mbr_dimensions
                                 : (time_series_size < 128 ? time_series_size : 128);
            if (trie_split_dimensions == 0) {
                trie_split_dimensions = index_segments < 32 ? index_segments : 32;
            }
            if (trie_split_dimensions < trie_bound_dimensions ||
                trie_split_dimensions > index_segments) {
                fprintf(stderr,
                        "error: trie split dimensions must be between record-bound dimensions (%d) and symbolic dimensions (%d).\n",
                        trie_bound_dimensions, index_segments);
                return EXIT_FAILURE;
            }
            if (function_type == 4 || function_type == 6) {
                if (n_coefficients == 0 || n_coefficients < index_segments) {
                    n_coefficients = index_segments;
                }
            }
        }


        if (!inmemory_flag) {
            sprintf(rm_command, "rm -rf %s", index_path);
            system(rm_command);
        }
        // check if n_segments size is at most timeseries_size
        if (index_segments > time_series_size) {
            fprintf(stderr, "ERROR: PAA segments may not be larger than timeseries-size!\n");
            return -1;
        }
        if (index_type == MESSI_INDEX_ISAX && time_series_size % n_segments != 0) {
            fprintf(stderr,
                    "WARNING: PAA ignores the final %d sample(s) because timeseries-size (%d) "
                    "is not divisible by n-segments (%d).\n",
                    time_series_size % n_segments, time_series_size, n_segments);
        }
        if (function_type == 4 || function_type == 6) {
            if (index_segments <= 0 || index_segments % 2 != 0) {
                fprintf(stderr, "ERROR: SFA/PISA n-segments must be a positive even number!\n");
                return -1;
            }
            if (n_coefficients != 0 &&
                (n_coefficients % 2 != 0 || n_coefficients < index_segments ||
                 n_coefficients > time_series_size)) {
                fprintf(stderr,
                        "ERROR: SFA/PISA coeff number must be even and between %d and %d!\n",
                        index_segments, time_series_size);
                return -1;
            }
            if (function_type == 4 && index_type == MESSI_INDEX_TRIE &&
                (n_coefficients == 0 || n_coefficients < index_segments ||
                 index_segments > time_series_size)) {
                fprintf(stderr,
                        "ERROR: trie SFA requires --sfa-n-coefficients (even, at least %d) and dimensions no greater than timeseries-size.\n",
                        index_segments);
                return -1;
            }
        }

        //get current time for logfiles
        time(&time_now);

        char time_str[20];
        time_now = time(NULL);
        strftime(time_str, 20, "%Y_%m_%d_%H:%M:%S", localtime(&time_now));

        //concatenate names for logfile directories
        char log_file_directory[FILENAME_LENGTH];
        char log_filename[FILENAME_LENGTH];
        char log_filename_tree[FILENAME_LENGTH];
        char log_filename_index[FILENAME_LENGTH];
        char log_filename_query[FILENAME_LENGTH];

        char default_log_root[FILENAME_LENGTH];
        const char *log_root = getenv("MESSI_LOG_ROOT");
        if (log_root == NULL || log_root[0] == '\0') {
            const char *home = getenv("HOME");
            if (home != NULL && home[0] != '\0' &&
                snprintf(default_log_root, sizeof(default_log_root), "%s/MESSI_logs", home) <
                    (int) sizeof(default_log_root)) {
                log_root = default_log_root;
            } else {
                log_root = "MESSI_logs";
            }
        }

        snprintf(log_file_directory, sizeof(log_file_directory), "%s", log_root);
        snprintf(log_filename, sizeof(log_filename), "%s/settings", log_root);
        snprintf(log_filename_tree, sizeof(log_filename_tree), "%s/tree", log_root);
        snprintf(log_filename_index, sizeof(log_filename_index), "%s/index", log_root);
        snprintf(log_filename_query, sizeof(log_filename_query), "%s/query", log_root);

        if (ensure_directory(log_filename) != 0 ||
            ensure_directory(log_filename_tree) != 0 ||
            ensure_directory(log_filename_index) != 0 ||
            ensure_directory(log_filename_query) != 0) {
            fprintf(stderr, "warning: cannot create MESSI log directories below %s: %s\n",
                    log_root, strerror(errno));
        }
        //concatenate actual file names
        strcat(log_filename, "/MESSI_SETTINGS_");
        strcat(log_filename, time_str);
        strcat(log_filename, ".csv");

        strcat(log_filename_tree, "/MESSI_TREE_");
        strcat(log_filename_tree, time_str);
        strcat(log_filename_tree, ".csv");

        strcat(log_filename_index, "/MESSI_INDEX_");
        strcat(log_filename_index, time_str);
        strcat(log_filename_index, ".csv");

        strcat(log_filename_query, "/MESSI_QUERY_");
        strcat(log_filename_query, time_str);
        strcat(log_filename_query, ".csv");

        strcat(index_directory, time_str);

        //create logfiles
        FILE *logfile;
        logfile = open_logfile_or_tmp(log_filename);

        FILE *logfile_tree;
        logfile_tree = open_logfile_or_tmp(log_filename_tree);

        FILE *logfile_index;
        logfile_index = open_logfile_or_tmp(log_filename_index);

        FILE *logfile_query;
        logfile_query = open_logfile_or_tmp(log_filename_query);
        if (logfile == NULL || logfile_tree == NULL || logfile_index == NULL || logfile_query == NULL) {
            fprintf(stderr, "error: could not create required log streams.\n");
            return EXIT_FAILURE;
        }

        SET_LOGFILE(logfile_query);

        fprintf(logfile,
                "MESSI SETTINGS\nFunction type,%d\nSIMD,%u\ntimeseries length,%d\npaa segments,%d\nisax-cardinality,%d\nleaf size,%d\nsample-size,%d\nsample type,%d\n",
                function_type, SIMD_flag, time_series_size, n_segments, sax_cardinality, leaf_size, sample_size,
                sample_type);

        if (!inmemory_flag && (function_type == 3 || function_type == 4 || function_type == 5 || function_type == 6)) {
            fprintf(stderr, "warning: function_type %d requires in-memory mode; enabling --inmemory.\n",
                    function_type);
            inmemory_flag = 1;
        }

        isax_index_settings *index_settings = isax_index_settings_init(index_path,         // INDEX DIRECTORY
                                                                       time_series_size,   // TIME SERIES SIZE
                                                                       index_segments,     // symbolic dimensions
                                                                       sax_cardinality,    // SAX CARDINALITY IN BITS
                                                                       leaf_size,          // LEAF SIZE
                                                                       min_leaf_size,      // MIN LEAF SIZE
                                                                       initial_lbl_size,   // INITIAL LEAF BUFFER SIZE
                                                                       flush_limit,        // FLUSH LIMIT
                                                                       initial_fbl_size,   // INITIAL FBL BUFFER SIZE
                                                                       total_loaded_leaves,// Leaves to load at each fetch
                                                                       tight_bound,        // Tightness of leaf bounds
                                                                       aggressive_check,   // aggressive check
                                                                       1,         // new index
                                                                       function_type,      //function_type
                                                                       inmemory_flag,      //inmemory_flag
                                                                       SIMD_flag,          //SIMD_flag
                                                                       sample_size,        //sample_size for MCB
                                                                       is_norm,            //input normalized for fft
                                                                       histogram_type,     //histogram type for binning
                                                                       sample_type,        //sampling type
                                                                       n_coefficients,      //coeff number
                                                                       index_type           //index layout
        );

        if (index_settings == NULL) {
            fprintf(stderr, "error: failed to initialize index settings.\n");
            return EXIT_FAILURE;
        }

        index_settings->node_split_criterion = node_split_criterion;
        index_settings->index_type = index_type;
        index_settings->trie_bound_dimensions = trie_bound_dimensions;
        index_settings->trie_split_dimensions = trie_split_dimensions;
        index_settings->trie_record_mbr_suffix_bound = trie_record_mbr_suffix_bound;
        index_settings->trie_fanout = trie_fanout;
        index_settings->trie_dynamic_alphabet = trie_dynamic_alphabet;
        index_settings->trie_min_bits = fanout_to_bits(trie_min_fanout);
        index_settings->trie_max_bits = fanout_to_bits(trie_max_fanout);
        index_settings->trie_alphabet_budget_bits = trie_alphabet_budget_bits;
        index_settings->sampling_seed = sampling_seed;
        index_settings->dynamic_root_split_variance =
                root_split_mode == MESSI_ROOT_SPLIT_VARIANCE;


        if (!inmemory_flag) {
            idx = isax_index_init(index_settings);
            print_settings(idx->settings, maxquerythread, trie_query_batch);
        } else {
            idx = isax_index_init_inmemory(index_settings);
            print_settings(idx->settings, maxquerythread, trie_query_batch);
        }

#ifdef CLUSTERED
        char s[255];
        sprintf(s, "rm -rf %s.*", dataset);
        system(s);
#endif

        COUNT_TOTAL_TIME_START

        //// #### COMMENTED OUT: BUGGY CODE!!! ####
        /*int sorted = 0;
        if(sorted == 1) {
        	isax_sorted_index_binary_file(dataset, dataset_size, idx);
        }
        else if(sorted == 2) {
        	isax_merge_sorted_index_binary_file(dataset, dataset_size, idx);
        }
        else {*/
        /// ########################################

        double query_wall_seconds = 0.0;
        if (inmemory_flag && index_type == MESSI_INDEX_TRIE && function_type == 3) {
            if (symbolic_trie_build(idx, dataset, dataset_size, filetype_int, apply_znorm) != SUCCESS) {
                fprintf(stderr, "error: trie construction failed.\n");
                return EXIT_FAILURE;
            }
            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);
            double query_wall_start = monotonic_seconds();
            if ((trie_query_batch ? symbolic_trie_query_file_batch : symbolic_trie_query_file)
                    (idx, queries, queries_size, filetype_int, apply_znorm, minimum_distance) != SUCCESS) {
                fprintf(stderr, "error: trie query processing failed.\n");
                return EXIT_FAILURE;
            }
            query_wall_seconds = monotonic_seconds() - query_wall_start;
            fprintf(stderr, ">>> query wall time: %.6f s\n", query_wall_seconds);

        //MESSI-SFA: in-memory flag set with function-type 4
        } else if (inmemory_flag && function_type == 4) {
            //initialize bins
            sfa_bins_init(idx);

            //set bins
            if (sfa_set_bins(idx, dataset, dataset_size, maxquerythread, filetype_int,
                             apply_znorm) != SUCCESS) {
                fprintf(stderr, "error: SFA bin preparation failed; aborting index creation.\n");
                return EXIT_FAILURE;
            }

            //build index
            if (index_type == MESSI_INDEX_TRIE) {
                if (symbolic_trie_build(idx, dataset, dataset_size, filetype_int, apply_znorm) != SUCCESS) {
                    fprintf(stderr, "error: trie construction failed.\n"); return EXIT_FAILURE;
                }
            } else index_creation_pRecBuf(dataset, dataset_size, filetype_int, apply_znorm, idx, dynamic_index);

            //calculate depth (for analysis logfile only)
            if (index_type == MESSI_INDEX_ISAX) calculate_average_depth(logfile_tree, idx);

            //save index building stats
            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);
            for (int i = 0; i < idx->settings->n_segments; i++) {
                memcpy(&idx->binsv[i * (idx->settings->sax_alphabet_cardinality - 1)], idx->bins[i],
                       sizeof(ts_type) * (idx->settings->sax_alphabet_cardinality - 1));
            }

            //perform queries
            // if (topk && k_size > 1) {
            //     isax_topk_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
            //                                             min_checked_leaves, k_size, filetype_int, apply_znorm,
            //                                             &exact_topk_MESSImq_inmemory);//MESSI topk
            // } else {
            double query_wall_start = monotonic_seconds();
            if (index_type == MESSI_INDEX_TRIE) {
                if ((trie_query_batch ? symbolic_trie_query_file_batch : symbolic_trie_query_file)(idx, queries, queries_size, filetype_int, apply_znorm, minimum_distance) != SUCCESS) return EXIT_FAILURE;
            } else isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                       filetype_int, apply_znorm, dynamic_index, &exact_search_MESSI);
            query_wall_seconds = monotonic_seconds() - query_wall_start;
            fprintf(stderr, ">>> query wall time: %.6f s\n", query_wall_seconds);

        } else if (inmemory_flag && function_type == 5) {
            //initialize bins
            spartan_bins_init(idx);

            //set bins
            if (spartan_set_bins(idx, dataset, dataset_size, maxquerythread, filetype_int,
                                 apply_znorm) != SUCCESS) {
                fprintf(stderr, "error: SPARTAN bin preparation failed; aborting index creation.\n");
                return EXIT_FAILURE;
            }

            //build index
            if (index_type == MESSI_INDEX_TRIE) {
                if (symbolic_trie_build(idx, dataset, dataset_size, filetype_int, apply_znorm) != SUCCESS) {
                    fprintf(stderr, "error: trie construction failed.\n"); return EXIT_FAILURE;
                }
            } else index_creation_pRecBuf(dataset, dataset_size, filetype_int, apply_znorm, idx, dynamic_index);

            //calculate depth (for analysis logfile only)
            if (index_type == MESSI_INDEX_ISAX) calculate_average_depth(logfile_tree, idx);

            //save index building stats
            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);
            for (int i = 0; i < idx->settings->n_segments; i++) {
                memcpy(&idx->binsv[i * (idx->settings->sax_alphabet_cardinality - 1)], idx->bins[i],
                       sizeof(ts_type) * (idx->settings->sax_alphabet_cardinality - 1));
            }

            //perform queries
            /*if (topk && k_size > 1) {
                isax_topk_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
                                                        min_checked_leaves, k_size, filetype_int, apply_znorm,
                                                        &exact_topk_MESSImq_inmemory);//MESSI topk
            } else {*/
            double query_wall_start = monotonic_seconds();
            if (index_type == MESSI_INDEX_TRIE) {
                if ((trie_query_batch ? symbolic_trie_query_file_batch : symbolic_trie_query_file)(idx, queries, queries_size, filetype_int, apply_znorm, minimum_distance) != SUCCESS) return EXIT_FAILURE;
            } else isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                       filetype_int, apply_znorm, dynamic_index, &exact_search_MESSI);
            query_wall_seconds = monotonic_seconds() - query_wall_start;
            fprintf(stderr, ">>> query wall time: %.6f s\n", query_wall_seconds);

        } else if (inmemory_flag && function_type == 6) {
            //initialize bins
            pisa_bins_init(idx);

            //set bins
            if (pisa_set_bins(idx, dataset, dataset_size, maxquerythread, filetype_int,
                              apply_znorm) != SUCCESS) {
                fprintf(stderr, "error: PISA bin preparation failed; aborting index creation.\n");
                return EXIT_FAILURE;
            }

            //build index
            if (index_type == MESSI_INDEX_TRIE) {
                if (symbolic_trie_build(idx, dataset, dataset_size, filetype_int, apply_znorm) != SUCCESS) {
                    fprintf(stderr, "error: trie construction failed.\n"); return EXIT_FAILURE;
                }
            } else index_creation_pRecBuf(dataset, dataset_size, filetype_int, apply_znorm, idx, dynamic_index);

            //calculate depth (for analysis logfile only)
            if (index_type == MESSI_INDEX_ISAX) calculate_average_depth(logfile_tree, idx);

            //save index building stats
            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);
            for (int i = 0; i < idx->settings->n_segments; i++) {
                memcpy(&idx->binsv[i * (idx->settings->sax_alphabet_cardinality - 1)], idx->bins[i],
                       sizeof(ts_type) * (idx->settings->sax_alphabet_cardinality - 1));
            }

            //perform queries
            /*if (topk && k_size > 1) {
                isax_topk_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
                                                        min_checked_leaves, k_size, filetype_int, apply_znorm,
                                                        &exact_topk_MESSImq_inmemory);//MESSI topk
            } else {*/
            double query_wall_start = monotonic_seconds();
            if (index_type == MESSI_INDEX_TRIE) {
                if ((trie_query_batch ? symbolic_trie_query_file_batch : symbolic_trie_query_file)(idx, queries, queries_size, filetype_int, apply_znorm, minimum_distance) != SUCCESS) return EXIT_FAILURE;
            } else isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                       filetype_int, apply_znorm, dynamic_index, &exact_search_MESSI);
            query_wall_seconds = monotonic_seconds() - query_wall_start;
            fprintf(stderr, ">>> query wall time: %.6f s\n", query_wall_seconds);

        } else if (inmemory_flag) {
            // MESSI: parallel in memory index creation 
            index_creation_pRecBuf(dataset, dataset_size, filetype_int, apply_znorm, idx, dynamic_index);

            calculate_average_depth(logfile_tree, idx);

            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);

            //INIT_STATS()
            if (knnlabel) {
                if (function_type == 1) {
                    isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                               min_checked_leaves, k_size, 2000, &exact_search_serial_topk_inmemory);
                } else if (function_type == 2) {
                    isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                               min_checked_leaves, k_size, 2000, &exact_topk_serial_ParIS_inmemory);
                } else if (function_type == 3) {
                    isax_knn_query_binary_file_traditional(queries, labelset, queries_size, idx, minimum_distance,
                                                           min_checked_leaves, k_size, 2000,
                                                           &exact_topk_MESSImq_inmemory);
                }
            } /*else if (topk && k_size > 1) {
                if (function_type == 3)
                    isax_topk_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
                                                            min_checked_leaves, k_size, filetype_int, apply_znorm,
                                                            &exact_topk_MESSImq_inmemory);//MESSI topk

            } */ else {
                if (function_type == 0) {
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_inmemory);
                    //isax_DTWquery_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,dtwwindowsize);
                } else if (function_type == 1) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search_parads_inmemory);
                } else if (function_type == 4) {

                } else if (function_type == 7) {
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_nb_inmemory);
                } else if (function_type == 2) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                           &exact_search_serial_ParIS_inmemory);
                } else if (function_type == 5) {
                    //bf=1;
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS2_inmemory);
                    //isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, filetype_int, filetype_int, apply_znorm, &exact_search_ParISnew_inmemory_hybrid_workstealing);
                } else if (function_type == 6) {
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParGISG_openmp_inmemory);
                }
                    //MESSI-mq: in-memory flag set with function-type 3
                else if (function_type == 3) {
                    double query_wall_start = monotonic_seconds();
                    isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                       filetype_int, apply_znorm, dynamic_index, &exact_search_MESSI);
                    query_wall_seconds = monotonic_seconds() - query_wall_start;
                    fprintf(stderr, ">>> query wall time: %.6f s\n", query_wall_seconds);

                    //isax_query_binary_file_batch(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_ParIS_nb_batch_inmemory);
                } else if (function_type == 8) {
                    //isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, filetype_int, filetype_int, znorm, &exact_search_ParISnew_inmemory);
                    //isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, filetype_int, znorm, &exact_search_ParISnew_inmemory_workstealing);
                } else if (function_type == 9) {
                    // isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves, filetype_int, znorm, &exact_search_serial_ParGIS_openmp_inmemory);
                } else if (function_type == 10) {

                    //index_generate_inmemory_pRecBuf(dataset, dataset_size, idx);
                    //COUNT_QUEUE_TIME_START
                    //flush_pRecBuf_inmemory((parallel_first_buffer_layer*) idx->fbl, idx);
                    //COUNT_QUEUE_TIME_END
                    //PRINT_STATS(0.0f)

                    //printf("this is the queue time: %f\n", total_queue_time);
                    //INIT_STATS()
                    //isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial_inmemory);
                }
            }

        } else {
            // ParIS/ParIS+ on disk index creation
            if (function_type == 0) {
                // ADS+ index creation
                isax_index_binary_file(dataset, dataset_size, idx);

                flush_fbl(idx->fbl, idx);
            } else if (function_type == 1) {
                isax_index_binary_file_m(dataset, dataset_size, idx,
                                         calculate_thread);// ParIS indexing program in parallel
            } else if (function_type == 2) {
                isax_index_binary_file_m_new(dataset, dataset_size, idx,
                                             calculate_thread);//ParIS+ indexing program in parallel
            }
        }

        SAVE_STATS_TOTAL(logfile_query, queries_size)
        if (queries_size > 0 && query_wall_seconds > 0.0) {
            const double avg_checked_nodes = (double) checked_nodes_all / queries_size;
            const double avg_lower_bounds = (double) LBDcalculationnumber_all / queries_size;
            const double avg_exact_distances = (double) RDcalculationnumber_all / queries_size;
            const double checked_node_percent = total_tree_nodes > 0
                                                ? 100.0 * avg_checked_nodes / total_tree_nodes
                                                : 0.0;
            const double lower_bound_percent = dataset_size > 0
                                               ? 100.0 * avg_lower_bounds / dataset_size
                                               : 0.0;
            const double exact_distance_percent = dataset_size > 0
                                                  ? 100.0 * avg_exact_distances / dataset_size
                                                  : 0.0;
            /* In a trie every series reaching a leaf receives exactly one
             * record lower-bound evaluation.  This lets the summary separate
             * pruning before leaf scans (node MBRs) from pruning by the
             * record bound itself, without adding counters to hot loops. */
            const double node_mbr_skipped = avg_lower_bounds < dataset_size
                                                ? (double) dataset_size - avg_lower_bounds
                                                : 0.0;
            const double record_bound_pruned = avg_exact_distances < avg_lower_bounds
                                                   ? avg_lower_bounds - avg_exact_distances
                                                   : 0.0;
            const double node_mbr_percent = dataset_size > 0
                                                ? 100.0 * node_mbr_skipped / dataset_size
                                                : 0.0;
            const double record_bound_candidate_percent = avg_lower_bounds > 0.0
                                                              ? 100.0 * record_bound_pruned / avg_lower_bounds
                                                              : 0.0;
            const double record_bound_index_percent = dataset_size > 0
                                                          ? 100.0 * record_bound_pruned / dataset_size
                                                          : 0.0;
            const double total_pruned_percent = dataset_size > 0
                                                   ? 100.0 * ((double) dataset_size - avg_exact_distances) /
                                                         dataset_size
                                                   : 0.0;
            char nodes[32], lower_bounds[32], exact_distances[32], index_nodes[32], indexed_series[32];
            char wall_time[32];
            format_compact_count(avg_checked_nodes, nodes, sizeof(nodes));
            format_compact_count(avg_lower_bounds, lower_bounds, sizeof(lower_bounds));
            format_compact_count(avg_exact_distances, exact_distances, sizeof(exact_distances));
            format_compact_count(total_tree_nodes, index_nodes, sizeof(index_nodes));
            format_compact_count(dataset_size, indexed_series, sizeof(indexed_series));
            if (query_wall_seconds < 1.0)
                snprintf(wall_time, sizeof(wall_time), "%.3f ms", 1000.0 * query_wall_seconds);
            else
                snprintf(wall_time, sizeof(wall_time), "%.3f s", query_wall_seconds);
            fprintf(stderr, "=== Query summary ===\n"
                   "  queries          : %d\n"
                   "  wall time        : %s (%.3f ms/query)\n"
                   "  checked nodes    : %s/query (%.2f%% of %s index nodes)\n"
                   "  lower bounds     : %s/query (%.2f%% of %s indexed series)\n"
                   "  exact distances  : %s/query (%.2f%% of %s indexed series)\n",
                   queries_size, wall_time, 1000.0 * query_wall_seconds / queries_size,
                   nodes, checked_node_percent, index_nodes,
                   lower_bounds, lower_bound_percent, indexed_series,
                   exact_distances, exact_distance_percent, indexed_series);
            if (index_type == MESSI_INDEX_TRIE) {
                const char *record_bound_name = idx->settings->trie_record_mbr_suffix_bound
                    ? "prefix + MBR suffix" : "symbolic record bound";
                fprintf(stderr, "  pruning breakdown:\n"
                       "    %-20s : %.2f%% indexed series skipped before record bounds\n"
                       "    %-20s : %.2f%% of reached records pruned (%.2f%% of indexed series)\n"
                       "    %-20s : %.2f%% indexed series pruned before exact distance\n",
                       "node MBRs", node_mbr_percent,
                       record_bound_name, record_bound_candidate_percent,
                       record_bound_index_percent,
                       "total", total_pruned_percent);
            }
            if (profile_query_phases) {
                fprintf(stderr, "  phase profile (accumulated worker ms/query):\n"
                       "    node MBR bounds  : %.3f\n"
                       "    record bounds    : %.3f\n"
                       "    exact distances  : %.3f\n",
                       total_mbr_dist_calc_time_all / (1000.0 * queries_size),
                       total_record_lb_dist_calc_time_all / (1000.0 * queries_size),
                       total_real_dist_calc_time_all / (1000.0 * queries_size));
                if (index_type == MESSI_INDEX_TRIE) {
                    fprintf(stderr, "    frontier traversal: %.3f\n"
                           "    queue locks/pops : %.3f\n"
                           "    candidate heap   : %.3f\n"
                           "    synchronization/wait: %.3f\n",
                           total_trie_frontier_time_all / (1000.0 * queries_size),
                           total_trie_queue_time_all / (1000.0 * queries_size),
                           total_trie_heap_time_all / (1000.0 * queries_size),
                           total_trie_sync_time_all / (1000.0 * queries_size));
                } else {
                    fprintf(stderr, "    traversal/queues : %.3f\n",
                           total_tree_pass_time_all / (1000.0 * queries_size));
                }
            }
        }

        /* Do not store index for now
        //save index and get size for analysis
        index_mRecBuf_write(idx);

        struct stat stat_index;
        struct stat stat_adaptive;

        get_index_size(idx, &stat_index, &stat_adaptive);

        fprintf(stderr, "\nindex size: %ld\n", (long int) stat_index.st_size);
        fprintf(logfile_index, "%ld\n", (long int) stat_index.st_size);
        */

        fclose(logfile);
        fclose(logfile_tree);
        fclose(logfile_query);
        fclose(logfile_index);

        float distance = 0;
        FLUSHES++;

        COUNT_TOTAL_TIME_END
        //PRINT_STATS(distance)

        COUNT_TOTAL_TIME_START
        if (complete_type == 1) {
            fprintf(stderr, ">>> Completing index.\n");
            complete_index(idx, dataset_size);
        } else if (complete_type == 2) {
            fprintf(stderr, ">>> Completing index.\n");
            complete_index_leafs(idx);
        }
        COUNT_TOTAL_TIME_END

        if (!inmemory_flag) {
            index_write(idx);
            isax_index_destroy(idx, NULL);
            PRINT_STATS(distance)

        } else {
            free(rawfile);
            rawfile = NULL;
            if (index_type == MESSI_INDEX_TRIE) symbolic_trie_destroy(idx);

            if (function_type == 4) {
                sfa_free_bins(idx);
                isax_index_pRecBuf_destroy(idx, NULL, maxquerythread);
            } else if (function_type == 5) {
                spartan_free_bins(idx);
                isax_index_pRecBuf_destroy(idx, NULL, maxquerythread);
            } else if (function_type == 6) {
                pisa_free_bins(idx);
                isax_index_pRecBuf_destroy(idx, NULL, maxquerythread);
            } else if (function_type == 3) {
                isax_index_pRecBuf_destroy(idx, NULL, maxquerythread);
            } else {
                //isax_index_destroy(idx, NULL);
                MESSI2_index_destroy(idx, NULL);
            }
        }

        fprintf(stderr, "\n");

    }

    return 0;
}

void INThandler(int sig) {
    static const char message[] = "\nInterrupted.\n";

    (void)sig;
    /* write(2) and _exit(2) are safe in a signal handler.  In particular, do
     * not prompt, acquire index locks, or attempt to save a partial index
     * while a worker may be modifying it. */
    (void)write(STDERR_FILENO, message, sizeof(message) - 1);
    _exit(130);
}
