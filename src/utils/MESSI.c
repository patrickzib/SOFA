//
//  MESSI.c
//  isaxlib
//
//  created: Jakob Brand 2022
//  changes: Gerrit Slomma 2025
//
#define _GNU_SOURCE

#define PRODUCT "----------------------------------------------\
\nThis is the Adaptive Leaf iSAX index.\n\
----------------------------------------------\n\n"
#ifdef VALUES
#include <values.h>
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <sched.h>
#if defined(__linux__)
#include <numa.h>
#endif
#if defined(__APPLE_CC__)
#include <sys/sysctl.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>

#include "ads/sax/sax.h"
#include "ads/sax/ts.h"
#include "ads/isax_visualize_index.h"
#include "ads/isax_file_loaders.h"
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
#include "include/ads/isax_file_loaders.h"

#include "../../config.h"
#include "../../globals.h"

isax_index *idx;

void INThandler(int);

int main(int argc, char **argv) {
    signal(SIGINT, INThandler);

#if defined(__x86_64__)
    if (!__builtin_cpu_supports("fma")) {
        printf("CPU does not support FMA, exiting.\n");

        return -1;
    }
#endif

#ifndef BENCHMARK
    printf("%s", PRODUCT);

#if VERBOSE_LEVEL == 0
    printf("Executing in silent mode. Please wait.\n");
#endif
#endif
    static char index_directory[FILENAME_MAX];
    int pathlen = snprintf(index_directory, strlen(getenv("HOME")) + 1, getenv("HOME"));
    snprintf(index_directory + pathlen, strlen("/index/") + 1, "/index/");

    static char data_directory[FILENAME_MAX];
    pathlen = snprintf(data_directory, strlen(getenv("HOME")) + 1, getenv("HOME"));
    snprintf(data_directory + pathlen, strlen("/data/data") + 1, "/data/data");

    static char query_directory[FILENAME_MAX];
    pathlen = snprintf(query_directory, strlen(getenv("HOME")) + 1, getenv("HOME"));
    snprintf(query_directory + pathlen, strlen("/data/query") + 1, "data/query");

    static char *dataset = data_directory;
    static char *queries = query_directory;
    static char *index_path = index_directory;
    static char *labelset = index_directory;
    static int64_t dataset_size = 6000000;  // testbench
    static int queries_size = 10;
    static int time_series_size = 256;
    static int paa_segments = 16;
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
    static int cpu_control_type = 81;
    static int threads_to_run = 8;
    static int node_to_run = -1;
    static char inmemory_flag = 0;
    static char SIMD_flag = 0;
    static char is_norm = 0;
    static int histogram_type = 1;
    static int sample_type = 1;
    static int coeff_number = 0;
    static int filetype_int = 0;
    static int apply_znorm = 0;

    int calculate_thread = 8;
    int function_type = 0;
    N_PQUEUE = 1;
    maxreadthread = 5;
    read_block_length = 20000;
    int k_size = 0;
    int64_t labelsize = 1;
    int topk = 0;
    int dtwwindowsize = 0;
    int sample_size = 1;

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
                {"sax-cardinality",     required_argument, 0,    'x'},
                {"paa-segments",        required_argument, 0,    'B'},
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
                {"coeff-number",        required_argument, 0,    'D'},
                {"filetype-int",        no_argument,       0,    'E'},
                {"apply-z-norm",        no_argument,       0,    'F'},
                {NULL,                  0,                 NULL, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;
        int c = getopt_long(argc, argv, "", long_options, &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
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
                paa_segments = atoi(optarg);

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
                if (strchr(optarg, '-')) {
                    node_to_run = atoi(strchr(optarg, '-') + 1);
                }

                threads_to_run = atoi(optarg);

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
                coeff_number = atoi(optarg);

                break;
            case 'E':
                filetype_int = 1;

                break;
            case 'F':
                apply_znorm = 1;

                break;
            case 'h':
#ifdef BENCHMARK
                printf("%s", PRODUCT);
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
                \t--paa-segments XX\t\tThe number of segments to divide the time series.\n\
                \t--leaf-size XX\t\t\tThe maximum size of each leaf\n\
                \t--min-leaf-size XX\t\tThe minimum size of each leaf\n\
                \t--initial-lbl-size XX\t\tThe initial lbl buffer size for each buffer.\n\
                \t--flush-limit XX\t\tThe limit of time series in memory at the same time\n\
                \t--initial-fbl-size XX\t\tThe initial fbl buffer size for each buffer.\n\
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
                \t\t\tMESSI-mq: 3\n\
                \t\t\tMESSI-SFA: 4\n\
                \t--SIMD\t\t\tSet for search with SIMD intrinsics\n\
                \t--sample-size\t\t\tSet sample size for MCB\n\
                \t--sample-type\t\t\tSet for sampling strategy\n\
                \t\t\tfirst-n-values sampling: 1\n\
                \t\t\tuniform sampling: 2\n\
                \t\t\trandom sampling: 3\n\
                \t--filetype-int\t\t\tSet if the input time series file is stored in int-type\n\
                \t--apply-z-norm\t\t\tApply z-normalization to the data\n\
                \t--is-norm\t\t\tSet for search with normalized input time series\n\
                \t--coeff-number\t\t\tSet number of coeff to choose highest-variance coeff (doubled for real & imag parts - must be between paa_segments/2 and timeseries-size/2)\n\
                \t--histogram-type\t\t\tSet for binning strategy\n\
                \t\t\tequi-depth splitting (default): 1\n\
                \t\t\tequi-width splitting: 2\n\
                \t--cpu-type\t\t\tSet for how many cores you want to used on what cpu numa-node\n\
                \t--help\n\n\
                \tCPU type code:\t\t\t2-1 : 2 core on node 1\n\
                \t\t\t\t\t2 : 2 cores in all nodes\n\
                \t\t\t\t\t40-0 : 40 cores on node 0\n");

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

#if defined(__linux__)
    if (numa_available() >= 0) {
        struct bitmask *bm = numa_bitmask_alloc(numa_num_task_cpus());
        numa_set_localalloc();

        if (node_to_run >= 0) {
            if (node_to_run > numa_max_node()) {
                node_to_run = numa_max_node();
            }

            numa_run_on_node(node_to_run);
            printf("run %d threads on node %d\n", threads_to_run, node_to_run);
        } else {
            printf("run %d threads on all nodes\n", threads_to_run);
            numa_run_on_node(-1);
        }
    } else {
#endif
        printf("run %d threads on all nodes\n", threads_to_run);
#if defined(__linux__)
    }
#endif

    calculate_thread = threads_to_run;
    maxquerythread = threads_to_run;

    if (use_index) {
        isax_index *idx = index_read(index_path);
        idx->settings->tight_bound = tight_bound;
        idx->settings->aggressive_check = aggressive_check;
        idx->settings->total_loaded_leaves = total_loaded_leaves;
        idx->settings->min_leaf_size = min_leaf_size;
        print_settings(idx->settings);

        char sanity_test = 0;

        if (sanity_test) {
            cache_sax_file(idx);
            isax_query_binary_file(queries, queries_size, idx,
                    minimum_distance, min_checked_leaves, &sanity_check_query);
        } else {
            if (serial_scan) {
                cache_sax_file(idx);
                if (knnlabel) {
                    if (function_type == 0) {
                        isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                min_checked_leaves, k_size, labelsize, &exact_topk_serial);  // ADS+ knn
                    } else if (function_type == 1) {
                        isax_knn_query_binary_file(queries, labelset, queries_size, idx, minimum_distance,
                                min_checked_leaves, k_size, labelsize, &exact_topk_serial_ParIS);  // ParIS knn
                    }
                } else if (topk) {
                    if (function_type == 0) {
                        isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, k_size, &exact_topk_serial);  // ADS+ topk
                    } else if (function_type == 1) {
                        isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, k_size, &exact_topk_serial_ParIS);  // ParIS topk
                    }
                } else {
                    if (function_type == 0) {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, &exact_search_serial);  // ADS+ similarity search
                    } else if (function_type == 1) {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, &exact_search_serial_ParIS);  // ParIS similarity search
                    } else if (function_type == 2) {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, &exact_search_serial_ParIS_nb);  // ParIS-nb similarity search
                    } else if (function_type == 3) {
                        // TODO(someon): implement serial scan
                        // isax_query_binary_file_para(queries, queries_size, idx, minimum_distance, min_checked_leaves, &exact_search_serial);
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
                        isax_topk_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, k_size, &exact_topk);
                    }
                } else {
                    if (function_type == 0) {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, &exact_search);
                    } else {
                        isax_query_binary_file(queries, queries_size, idx, minimum_distance,
                                min_checked_leaves, &exact_search_m);
                    }
                }
            }
        }

        PRINT_STATS(0.00f)

        flush_all_leaf_buffers(idx, TMP_ONLY_CLEAN);
        index_write(idx);

        isax_index_destroy(idx, NULL);
    } else {
        char rm_command[FILENAME_MAX];

        if (!inmemory_flag) {
            snprintf(rm_command, FILENAME_MAX, "rm -rf %s", index_path);
            system(rm_command);
        }

        // check if paa_segments size is at most timeseries_size
        if (paa_segments > time_series_size) {
            fprintf(stderr, "ERROR: PAA segments may not be larger than timeseries-size!\n");

            return -1;
        }

        // check is coeff_number is between paa_segments/2 and timeseries_size/2
        if (coeff_number != 0 && (coeff_number < paa_segments / 2 || coeff_number > time_series_size / 2)) {
            if (coeff_number < paa_segments || coeff_number > time_series_size) {
                fprintf(stderr, "ERROR: coeff number must be between %d and %d!\n",
                        paa_segments, time_series_size);

                return -1;
            } else if (coeff_number % 2 != 0) {
                fprintf(stderr, "ERROR: coeff number must be divisible by 2!\n");

                return -1;
            }
        }

        char cpu_model[100];
        char cpu_cores[8];
#if defined(__linux__)
        FILE *fp = fopen("/proc/cpuinfo", "r");
        char buf[BUFSIZ];
        bool found_model = false;
        bool found_cores = false;

        while (fgets(buf, sizeof(buf), fp) != NULL) {
            if (strncmp(buf, "model name", strlen("model name")) == 0 && found_model == false) {
                char *pos = strstr(buf, ": ");

                if (pos != NULL) {
                    memcpy(cpu_model, pos + 2, strlen(pos) - 3);
                    found_model = true;
                }
            }

            if (strncmp(buf, "cpu cores", strlen("cpu cores")) == 0 && found_cores == false) {
                char *pos = strstr(buf, ": ");

                if (pos != NULL) {
                    memcpy(cpu_cores, pos + 2, strlen(pos) - 3);
                    found_cores = true;
                }
            }
        }

        printf("model: %s\n", cpu_model);
        printf("cores: %s\n", cpu_cores);
        fclose(fp);
#elif defined(__APPLE_CC__)
        int mib[2];
        size_t len;
        int cpu_cores_i;

        mib[0] = CTL_HW;
        mib[1] = HW_MODEL;

        sysctl(mib, 2, NULL, &len, NULL, 0);
        sysctl(mib, 2, cpu_model, &len, NULL, 0);

        mib[1] = HW_NCPU;

        sysctl(mib, 2, NULL, &len, NULL, 0);

        sysctl(mib, 2, &cpu_cores_i, &len, NULL, 0);
        sprintf(cpu_cores, "%d", cpu_cores_i);
#endif

        // get current time for logfiles
        time(&time_now);

        char time_str[20];
        time_now = time(NULL);
        strftime(time_str, 20, "%Y_%m_%d_%H:%M:%S", localtime(&time_now));

        // concatenate names for logfile directories
        char log_file_directory[FILENAME_MAX];
        char log_filename[FILENAME_MAX];
        char log_filename_tree[FILENAME_MAX];
        char log_filename_index[FILENAME_MAX];
        char log_filename_query[FILENAME_MAX];

        strcat(strcpy(log_file_directory, getenv("HOME")), "/MESSI_logs");
        strcat(strcpy(log_filename, getenv("HOME")), "/MESSI_logs/settings");
        strcat(strcpy(log_filename_tree, getenv("HOME")), "/MESSI_logs/tree");
        strcat(strcpy(log_filename_index, getenv("HOME")), "/MESSI_logs/index");
        strcat(strcpy(log_filename_query, getenv("HOME")), "/MESSI_logs/query");

        // check if logfile directories exist, create them if neccessary
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

        // concatenate actual file names
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

        // create logfiles
        FILE *logfile;
        logfile = fopen(log_filename, "w");

        FILE *logfile_tree;
        logfile_tree = fopen(log_filename_tree, "w");

        FILE *logfile_index;
        logfile_index = fopen(log_filename_index, "w");

        FILE *logfile_query;
        logfile_query = fopen(log_filename_query, "w");

        SET_LOGFILE(logfile_query);

        fprintf(logfile,
                "MESSI SETTINGS\nFunction type,%d\nSIMD,%u\ntimeseries length,%d\npaa segments,%d\nisax-cardinality,%d\nleaf size,%d\nsample-size,%d\nsample type,%d\nthreads run,%d\ncpu model,%s\ncpu cores,%s\n",
                function_type, SIMD_flag, time_series_size, paa_segments, sax_cardinality,
                leaf_size, sample_size, sample_type, threads_to_run, cpu_model, cpu_cores);

        isax_index_settings *index_settings = isax_index_settings_init(
                index_path,           // INDEX DIRECTORY
                time_series_size,     // TIME SERIES SIZE
                paa_segments,         // PAA SEGMENTS
                sax_cardinality,      // SAX CARDINALITY IN BITS
                leaf_size,            // LEAF SIZE
                min_leaf_size,        // MIN LEAF SIZE
                initial_lbl_size,     // INITIAL LEAF BUFFER SIZE
                flush_limit,          // FLUSH LIMIT
                initial_fbl_size,     // INITIAL FBL BUFFER SIZE
                total_loaded_leaves,  // Leaves to load at each fetch
                tight_bound,          // Tightness of leaf bounds
                aggressive_check,     // aggressive check
                1,                    // new index
                function_type,        // function_type
                inmemory_flag,        // inmemory_flag
                SIMD_flag,            // SIMD_flag
                sample_size,          // sample_size for MCB
                is_norm,              // input normalized for fft
                histogram_type,       // histogram type for binning
                sample_type,          // sampling type
                coeff_number);        // coeff number

        if (!inmemory_flag) {
            idx = isax_index_init(index_settings);
            print_settings(idx->settings);
        } else {
            idx = isax_index_init_inmemory(index_settings);
            print_settings(idx->settings);
        }

#ifdef CLUSTERED
        char s[255];
        snprintf(s, sizeof(s) + 1, "rm -rf %s.*", dataset);
        system(s);
#endif

#if defined(__x86_64__)
        if (__builtin_cpu_supports("fma")) {
            printf("FMA");
        }

        if (__builtin_cpu_supports("avx512f")) {
            printf("; AVX512F");
        }

        if (__builtin_cpu_supports("avx512vl")) {
            printf("; AVX512VL");
        }

        if (__builtin_cpu_supports("avx512dq")) {
            printf("; AVX512DQ");
        }

        printf(" detected\n");
#endif

        COUNT_TOTAL_TIME_START

        // MESSI-SFA: in-memory flag set with function-type 4
        if (inmemory_flag && function_type == 4) {
            // initialize bins
            sfa_bins_init(idx);

            // set bins
            sfa_set_bins(idx, dataset, dataset_size, maxquerythread, filetype_int, apply_znorm);

            // build index
            index_creation_pRecBuf(dataset, dataset_size, filetype_int, apply_znorm, idx);

            // calculate depth (for analysis logfile only)
            calculate_average_depth(logfile_tree, idx);

            // save index building stats
            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);

            for (int i = 0; i < paa_segments; i++) {
                memcpy(&idx->binsv[i * (idx->settings->sax_alphabet_cardinality - 1)], idx->bins[i],
                       sizeof(ts_type) * (idx->settings->sax_alphabet_cardinality - 1));
            }

            // perform queries
            if (topk && k_size > 1) {
                isax_topk_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
                        min_checked_leaves, k_size, filetype_int, apply_znorm,
                        &exact_topk_MESSImq_inmemory);  // MESSI topk
            } else {
                isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
                        min_checked_leaves, filetype_int, apply_znorm, &exact_search_MESSI);
            }
        } else if (inmemory_flag) {
            // MESSI: parallel in memory index creation
            index_creation_pRecBuf(dataset, dataset_size, filetype_int, apply_znorm, idx);

            calculate_average_depth(logfile_tree, idx);

            INIT_INDEX_STATS_FILE(logfile_index);
            INIT_SAVE_FILE(logfile_query);

            // INIT_STATS()
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
            } else if (topk && k_size > 1) {
                if (function_type == 3) {
                    isax_topk_query_binary_file_traditional(queries, queries_size, idx, minimum_distance,
                            min_checked_leaves, k_size, filetype_int, apply_znorm,
                            &exact_topk_MESSImq_inmemory);  // MESSI topk
                }
            } else {
                if (function_type == 1) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                            &exact_search_parads_inmemory);
                } else if (function_type == 2) {
                    isax_query_binary_file(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                            &exact_search_serial_ParIS_inmemory);
                } else if (function_type == 3) {
                    // MESSI-mq: in-memory flag set with function-type 3
                    isax_query_binary_file_traditional(queries, queries_size, idx, minimum_distance, min_checked_leaves,
                            filetype_int, apply_znorm, &exact_search_MESSI);
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
                            calculate_thread);   // ParIS indexing program in parallel
            } else if (function_type == 2) {
                isax_index_binary_file_m_new(dataset, dataset_size, idx,
                            calculate_thread);   // ParIS+ indexing program in parallel
            }
        }

        // save querying stats
        SAVE_STATS_TOTAL(logfile_query, queries_size)
        PRINT_STATS(0.00f)

        fclose(logfile);
        fclose(logfile_tree);
        fclose(logfile_query);
        fclose(logfile_index);

        float distance = 0;
        FLUSHES++;

        COUNT_TOTAL_TIME_END

        COUNT_TOTAL_TIME_START

        if (complete_type == 1) {
            fprintf(stderr, ">>> Completing index.\n");
            complete_index(idx, dataset_size);
        } else if (complete_type == 2) {
            fprintf(stderr, ">>> Completing index.\n");
            complete_index_leafs(idx);
        }

        COUNT_TOTAL_TIME_END

//        index_mRecBuf_write(idx);

//        for (int i = 0; i < paa_segments; i++) {
//            printf("[%2d] ", i);
//
//            for (int j = 0; j < idx->settings->sax_alphabet_cardinality - 1; j++) {
//                printf("%f ", idx->bins[i][j]);
//            }
//
//            printf("\n");
//        }

        if (!inmemory_flag) {
            index_write(idx);
            isax_index_destroy(idx, NULL);
            PRINT_STATS(distance)
        } else {
            free(rawfile);

            if (function_type == 4) {
                sfa_free_bins(idx);
                isax_index_pRecBuf_destroy(idx, NULL, maxquerythread);
            } else if (function_type == 3) {
                isax_index_pRecBuf_destroy(idx, NULL, maxquerythread);
            } else {
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
    } else {
        signal(SIGINT, INThandler);
    }

    getchar();  // Get new line character
}
