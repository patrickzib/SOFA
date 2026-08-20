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

void INThandler(int);

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

static double monotonic_seconds(void) {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (double) now.tv_sec + (double) now.tv_nsec / 1000000000.0;
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
    static int cpu_control_type = 0;
    static int requested_threads = 0;
    static int requested_numa_nodes = -1;
    static int threads_option_set = 0;
    static int numa_option_set = 0;
    static char inmemory_flag = 0;
    static char SIMD_flag = 0;
    static char is_norm = 0;
    static int histogram_type = 1;
    static int sample_type = 1;
    static int n_coefficients = 0;
    static int filetype_int = 0;
    static int apply_znorm = 0;
    static int dynamic_index = 1;
    static int dynamic_root_split_variance = 0;
    static int dynamic_root_split_uniform_set = 0;
    static int node_split_criterion = 1;
    static messi_index_type index_type = MESSI_INDEX_ISAX;
    static int symbolic_trie_dimensions = 0;
    /* Trie queries are independent and batch scheduling preserves each
     * query's private best-so-far, unlike speculative subtree parallelism. */
    static int trie_query_batch = 1;

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
                {"cpu-type",            required_argument, 0,    'w'},
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
                {"SIMD",                no_argument,       0,    '7'},
                {"sample-size",         required_argument, 0,    '8'},
                {"is-norm",             no_argument,       0,    '9'},
                {"histogram-type",      required_argument, 0,    'A'},
                {"sample-type",         required_argument, 0,    'C'},
                {"sfa-n-coefficients",  required_argument, 0,    'D'},
                {"filetype-int",        no_argument,       0,    'E'},
                {"apply-z-norm",        no_argument,       0,    'F'},
                {"node-split-criterion",required_argument, 0,   'G'},
                {"dynamic-index",       required_argument, 0,   'H'},
                {"dynamic-root-split-uniform", required_argument, 0, 'K'},
                {"dynamic-root-split-variance", no_argument,     0, 'L'},
                {"index-type", required_argument, 0, 1000},
                {"symbolic-trie-dimensions", required_argument, 0, 1001},
                {"trie-query-parallel", no_argument, 0, 1002},
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
                symbolic_trie_dimensions = atoi(optarg);
                break;
            case 1002:
                trie_query_batch = 0;
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
            case 'w':
                cpu_control_type = atoi(optarg);
                break;
            case 'I':
                threads_option_set = 1;
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
                numa_option_set = 1;
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
            case 'H':
                dynamic_index = atoi(optarg);
                dynamic_root_split_uniform_set = 1;
                break;
            case 'K':
                dynamic_index = atoi(optarg);
                dynamic_root_split_uniform_set = 1;
                break;
            case 'L':
                dynamic_root_split_variance = 1;
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
                \t--n-segments XX\t\tThe number of segments to divide the time series.\n\
                \t--leaf-size XX\t\t\tThe maximum size of each leaf\n\
                \t--min-leaf-size XX\t\tThe minimum size of each leaf\n\
                \t--initial-lbl-size XX\t\tThe initial lbl buffer size for each buffer.\n\
                \t--flush-limit XX\t\tThe limit of time series in memory at the same time\n\
                \t--initial-fbl-size XX\t\tThe initial fbl buffer size for each buffer.\n\
                \t--node-split-criterion XX\tSelect split decision (1=informed default, 2=simple, 3=maxvar, 4=maxbin)\n\
                \t--dynamic-index XX\t\tSet dynamic-index (kn) for root mask grouping\n\
                \t--dynamic-root-split-uniform XX\tUniform root split (alias for --dynamic-index)\n\
                \t--dynamic-root-split-variance\tUse all root-key bits, assigned by transform variance\n\
                \t--index-type isax|trie\tIndex layout (default: isax)\n\
                \t--symbolic-trie-dimensions N\tTrie symbolic dimensions, 16--64 (default: 64)\n\
                \t--trie-query-parallel\tParallelize each trie query across subtrees (default: batch queries)\n\
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
                \t--SIMD\t\t\tSet for search with SIMD intrinsics\n\
                \t--sample-size\t\t\tSet sample size for MCB\n\
                \t--sample-type\t\t\tSet for sampling strategy\n\
                \t\t\tfirst-n-values sampling: 1\n\
                \t\t\tuniform sampling: 2\n\
                \t\t\trandom sampling: 3\n\
                \t--filetype-int\t\t\tSet if the input time series file is stored in int-type\n\
                \t--apply-z-norm\t\t\tApply z-normalization to the data\n\
                \t--is-norm\t\t\tSet for search with normalized input time series\n\
                \t--sfa-n-coefficients\t\t\tSet number of coeff to choose highest-variance coeff (doubled for real & imag parts - must be between n_segments/2 and timeseries-size/2)\n\
                \t--histogram-type\t\t\tSet for binning strategy\n\
                \t\t\tequi-depth splitting (default): 1\n\
                \t\t\tequi-width splitting: 2\n\
                \t--threads N|auto\t\tWorker threads (default: all CPUs available to this process)\n\
                \t--numa auto|none|N\tCPU affinity: all available NUMA nodes, disabled, or first N nodes\n\
                \t--cpu-type\t\t\tDeprecated compatibility option; use --threads and --numa\n\
                \t--help\n\n\
                \tLegacy CPU type code:\t\tValues ending in 1 or 2 mean cores * 10 + NUMA nodes\n\
                \t\t\t\t\t22 : 2 core in 2 CPUs\n\
                \t\t\t\t\t41 : 4 core in 1 CPU\n\
                \t\t\t\t\t42 : 4 core in 2 CPUs\n\
                \t\t\t\t\t61 : 6 core in 1 CPU\n\
                \t\t\t\t\t62 : 6 core in 2 CPUs\n\
                \t\t\t\t\t81 : 8 core in 1 CPU\n\
                \t\t\t\t\t82 : 8 core in 2 CPUs\n\
                \t\t\t\t\t101: 10 core in 1 CPU\n\
                \t\t\t\t\t102: 10 core in 2 CPUs\n\
                \t\t\t\t\t121: 12 core in 1 CPU\n\
                \t\t\t\t\t122: 12 core in 2 CPUs\n\
                \t\t\t\t\t181: 18 core in 1 CPU\n\
                \t\t\t\t\t182: 18 core in 2 CPUs\n\
                \t\t\t\t\t242: 24 core in 2 CPUs\n\
                ");
                return 0;
                break;
            case 'a':
                use_index = 1;
                break;
            case 'v':
                inmemory_flag = 1;
                break;
            case '7':
                SIMD_flag = 1;
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

    /* --cpu-type is retained for existing invocations.  Its historical
     * encoding is cores * 10 + NUMA nodes, so 81 means 8 workers on one
     * node--not 81 workers.  --threads/--numa take precedence when supplied. */
    if (dynamic_root_split_variance && dynamic_root_split_uniform_set) {
        fprintf(stderr, "error: --dynamic-root-split-variance cannot be combined with a uniform root split.\n");
        return EXIT_FAILURE;
    }
    if (index_type == MESSI_INDEX_TRIE) {
        if (function_type < 4 || function_type > 6) {
            fprintf(stderr, "error: --index-type trie supports function types 4 (SFA), 5 (SPARTAN), and 6 (PISA).\n");
            return EXIT_FAILURE;
        }
        if (dynamic_root_split_uniform_set || dynamic_root_split_variance) {
            fprintf(stderr, "error: dynamic root split options are iSAX-only.\n");
            return EXIT_FAILURE;
        }
        if (symbolic_trie_dimensions == 0) symbolic_trie_dimensions = 64;
        if (symbolic_trie_dimensions < 16 || symbolic_trie_dimensions > 64) {
            fprintf(stderr, "error: symbolic-trie-dimensions must be between 16 and 64.\n");
            return EXIT_FAILURE;
        }
    }
    if (dynamic_index < 1 || dynamic_index > sax_cardinality) {
        fprintf(stderr, "error: uniform root split must be between 1 and sax-cardinality (%d).\n",
                sax_cardinality);
        return EXIT_FAILURE;
    }
    if (dynamic_root_split_variance &&
        function_type != 4 && function_type != 5 && function_type != 6) {
        fprintf(stderr, "error: --dynamic-root-split-variance is supported by SFA (4), SPARTAN (5), and PISA (6).\n");
        return EXIT_FAILURE;
    }

    if (cpu_control_type > 0 && !threads_option_set) {
        if (cpu_control_type >= 10 &&
            (cpu_control_type % 10 == 1 || cpu_control_type % 10 == 2)) {
            requested_threads = cpu_control_type / 10;
            if (!numa_option_set) {
                requested_numa_nodes = cpu_control_type % 10;
            }
        } else {
            requested_threads = cpu_control_type;
        }
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
    if (use_index) {
        isax_index *idx = index_read(index_path);
        idx->settings->tight_bound = tight_bound;
        idx->settings->aggressive_check = aggressive_check;
        idx->settings->total_loaded_leaves = total_loaded_leaves;
        idx->settings->min_leaf_size = min_leaf_size;
        print_settings(idx->settings);
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
        int index_segments = index_type == MESSI_INDEX_TRIE ? symbolic_trie_dimensions : n_segments;


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
                 n_coefficients > time_series_size / 2)) {
                fprintf(stderr,
                        "ERROR: SFA/PISA coeff number must be even and between %d and %d!\n",
                        index_segments, time_series_size / 2);
                return -1;
            }
            if (index_type == MESSI_INDEX_TRIE &&
                (n_coefficients == 0 || n_coefficients < index_segments ||
                 index_segments > time_series_size / 2)) {
                fprintf(stderr,
                        "ERROR: trie SFA/PISA requires --sfa-n-coefficients (even, at least %d) and dimensions no greater than timeseries-size/2.\n",
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

        strcat(strcpy(log_file_directory, getenv("HOME")), "/MESSI_logs");
        strcat(strcpy(log_filename, getenv("HOME")), "/MESSI_logs/settings");
        strcat(strcpy(log_filename_tree, getenv("HOME")), "/MESSI_logs/tree");
        strcat(strcpy(log_filename_index, getenv("HOME")), "/MESSI_logs/index");
        strcat(strcpy(log_filename_query, getenv("HOME")), "/MESSI_logs/query");

        //check if logfile directories exist, create them if neccessary
        struct stat st = {0};

        if (stat(log_file_directory, &st) == -1) {
            mkdir(log_file_directory, 0777);
        }
        if (stat(log_filename, &st) == -1) {
            mkdir(log_filename, 0777);
        }
        if (stat(log_filename_tree, &st) == -1) {
            mkdir(log_filename_tree, 0777);
        }
        if (stat(log_filename_index, &st) == -1) {
            mkdir(log_filename_index, 0777);
        }
        if (stat(log_filename_query, &st) == -1) {
            mkdir(log_filename_query, 0777);
        }
/*
        mkdir(log_filename, 0777);
        mkdir(log_filename_tree, 0777);
        mkdir(log_filename_index, 0777);
        mkdir(log_filename_query, 0777);
*/
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
                                                                       n_coefficients       //coeff number
        );

        if (index_settings == NULL) {
            fprintf(stderr, "error: failed to initialize index settings.\n");
            return EXIT_FAILURE;
        }

        index_settings->node_split_criterion = node_split_criterion;
        index_settings->index_type = index_type;
        index_settings->symbolic_trie_dimensions = symbolic_trie_dimensions;
        index_settings->dynamic_root_split_variance = dynamic_root_split_variance;


        if (!inmemory_flag) {
            idx = isax_index_init(index_settings);
            print_settings(idx->settings);
        } else {
            idx = isax_index_init_inmemory(index_settings);
            print_settings(idx->settings);
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

        //MESSI-SFA: in-memory flag set with function-type 4
        if (inmemory_flag && function_type == 4) {
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
            fprintf(stderr, ">>> query wall time: %.6f s\n", monotonic_seconds() - query_wall_start);

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
            fprintf(stderr, ">>> query wall time: %.6f s\n", monotonic_seconds() - query_wall_start);

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
            fprintf(stderr, ">>> query wall time: %.6f s\n", monotonic_seconds() - query_wall_start);

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
                    isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                                                       filetype_int, apply_znorm, dynamic_index, &exact_search_MESSI);

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

        //save querying stats
        if (inmemory_flag && (function_type == 4 || function_type == 5 || function_type == 6) &&
            queries_size > 0) {
            fprintf(stderr,
                    ">>> distance calculations: lower-bound=%lu (%.2f/query), exact=%lu (%.2f/query)\n",
                    LBDcalculationnumber_all,
                    (double) LBDcalculationnumber_all / queries_size,
                    RDcalculationnumber_all,
                    (double) RDcalculationnumber_all / queries_size);
        }
        SAVE_STATS_TOTAL(logfile_query, queries_size)
        PRINT_STATS(0.00f)

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

        printf("\n");

    }

    return 0;
}

void INThandler(int sig) {
    char c;

    signal(sig, SIG_IGN);
    fprintf(stderr, "Do you really want to quit? [y/n] ");
    c = getchar();
    if (c == 'y' || c == 'Y') {
        c = getchar();
        fprintf(stderr, "Do you want to save the index? [y/n] ");
        c = getchar();
        if (c == 'y' || c == 'Y') {
            flush_fbl(idx->fbl, idx);
            index_write(idx);
        }
        exit(0);
    } else
        signal(SIGINT, INThandler);
    getchar(); // Get new line character
}
