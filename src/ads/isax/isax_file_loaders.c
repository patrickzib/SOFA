//
//  isax_file_loaders.c
//  isax
//
//  Created by Kostas Zoumpatianos and Botao Peng, March 2020
//


#include "config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <pthread.h>
#include <semaphore.h>
#include <stdbool.h>
#include <float.h>
#include <unistd.h>
#include <math.h>

#include "ads/calc_utils.h"
#include "ads/isax_node.h"
#include "ads/isax_index.h"
#include "ads/isax_query_engine.h"
#include "ads/isax_node_record.h"
#include "ads/isax_file_loaders.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/inmemory_index_engine.h"
#include "ads/inmemory_query_engine.h"
#include "ads/sax/ts.h"
#include "ads/sfa/dft.h"
#include "ads/spartan/spartan.h"
#include "ads/pisa/pisa.h"

void isax_query_binary_file(const char *ifilename, int q_num, isax_index *index,
                            float minimum_distance, int min_checked_leaves,
                            query_result (*search_function)(ts_type *, ts_type *, isax_index *, float, int)) {
    fprintf(stderr, ">>> Performing queries in file: %s\n", ifilename);
    FILE *ifile;
    ifile = fopen(ifilename, "rb");
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }

    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);
    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    int q_loaded = 0;
    ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    ts_type *paa = malloc(sizeof(ts_type) * index->settings->n_segments);
    fftw_workspace fftw = {0};
    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        fftw_workspace_init(&fftw, index->settings->timeseries_size);
    }


    while (q_loaded < q_num) {
        COUNT_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_INPUT_TIME_END

        //Parse TS and make FFT representation
        if (index->settings->function_type == 4) {
            memcpy(fftw.ts, ts, sizeof(ts_type) * index->settings->timeseries_size);
            int use_best = index->settings->n_coefficients != 0;
            if (fft_from_ts(index, index->settings->n_segments, use_best, &fftw) != SUCCESS) {
                fprintf(stderr, "error: failed to calculate SFA query coefficients.\n");
                break;
            }

            memcpy(paa, fftw.transform, sizeof(ts_type) * index->settings->n_segments);
        } else if (index->settings->function_type == 5) {
            pca_from_ts(index, ts, paa);
        } else if (index->settings->function_type == 6) {
            pisa_pca_from_ts(index, ts, paa, &fftw);
        } else {
            paa_from_ts(ts, paa, index->settings);
        }

        COUNT_TOTAL_TIME_START
        query_result result = search_function(ts, paa, index, minimum_distance, min_checked_leaves);
        COUNT_TOTAL_TIME_END
        if (SHOULD_REPORT_QUERY(q_loaded, q_num)) PRINT_STATS(result.distance)

        fflush(stderr);
#if VERBOSE_LEVEL >= 1
        printf("[%p]: Distance: %lf\n", result.node, result.distance);
#endif

        q_loaded++;
    }
    free(paa);
    free(ts);
    fclose(ifile);
    fprintf(stderr, ">>> Finished querying.\n");
    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        fftw_workspace_destroy(&fftw);
    }

}

void isax_query_binary_file_traditional(
        const char *ifilename, int q_num, isax_index *index,
        float minimum_distance, int min_checked_leaves, int filetype_int,
        int apply_znorm, int kn,
        query_result (*search_function)(ts_type *, ts_type *, isax_index *, node_list *, float, int, int)) {

    fprintf(stderr, ">>> Performing queries in file: %s\n", ifilename);

    FILE *ifile;
    ifile = fopen(ifilename, "rb");
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }

    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);
    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    int q_loaded = 0;

    node_list nodelist;
    nodelist.nlist = malloc(sizeof(isax_node *) * pow(2, index->settings->n_segments));
    nodelist.node_amount = 0;
    isax_node *current_root_node = index->first_node;
    while (1) {
        if (current_root_node != NULL) {
            nodelist.nlist[nodelist.node_amount] = current_root_node;
            current_root_node = current_root_node->next;
            nodelist.node_amount++;
        } else {
            break;
        }

    }

    fprintf(stderr, ">>> node_amount is %d\n", nodelist.node_amount);

    fftw_workspace fftw = {0};

    file_type *ts_int32;
    if (filetype_int) {
        ts_int32 = malloc(sizeof(file_type) * index->settings->timeseries_size);
    }
    ts_type *paa = malloc(sizeof(ts_type) * index->settings->n_segments);
    ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    unsigned long ts_length = index->settings->timeseries_size;

    // create fftw plan
    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        fftw_workspace_init(&fftw, ts_length);
    }

    while (q_loaded < q_num) {
        COUNT_INPUT_TIME_START

        if (filetype_int) {
            fread(ts_int32, sizeof(file_type), ts_length, ifile);
            for (int i = 0; i < ts_length; ++i) {
                ts[i] = (ts_type) ts_int32[i];
            }

        } else {
            fread(ts, sizeof(ts_type), ts_length, ifile);
        }

        // apply z-normalization
        if (apply_znorm) {
            znorm(ts, ts_length);
        }

        COUNT_INPUT_TIME_END
        int print_query_row = SHOULD_REPORT_QUERY(q_loaded, q_num);
        if (print_query_row) {
            PRINT_STATS_HEADER();
            fprintf(stderr, "%3d: ", q_loaded);
        }

        COUNT_QUERYING_TIME_START
        COUNT_INIT_TIME_START

        if (index->settings->function_type == 4) {
            //SFA: parse ts and make fft representation
            memcpy(fftw.ts, ts, sizeof(ts_type) * ts_length);

            int use_best = index->settings->n_coefficients != 0;
            if (fft_from_ts(index, index->settings->n_segments, use_best, &fftw) != SUCCESS) {
                fprintf(stderr, "error: failed to calculate SFA query coefficients.\n");
                break;
            }

            memcpy(paa, fftw.transform, sizeof(ts_type) * index->settings->n_segments);
        } else if (index->settings->function_type == 5) {
            pca_from_ts(index, ts, paa);
        } else if (index->settings->function_type == 6) {
            pisa_pca_from_ts(index, ts, paa, &fftw);
        } else {
            // Parse ts and make PAA representation
            paa_from_ts(ts, paa, index->settings);
        }

        COUNT_TOTAL_TIME_START
        query_result result = search_function(ts, paa, index, &nodelist, minimum_distance, min_checked_leaves, kn);
        COUNT_TOTAL_TIME_END
        COUNT_QUERYING_TIME_END

        if (print_query_row) PRINT_STATS(result.distance)
        SAVE_STATS(result.distance)

        fflush(stderr);
#if VERBOSE_LEVEL >= 1
        printf("[%p]: Distance: %lf\n", result.node, result.distance);
#endif
        q_loaded++;
    }

    if (index->settings->function_type == 4 || index->settings->function_type == 6) {
        fftw_workspace_destroy(&fftw);
    }

    free(nodelist.nlist);
    free(paa);
    free(ts);

    if (filetype_int) {
        free(ts_int32);
    }
    fclose(ifile);

    fprintf(stderr, ">>> Finished querying.\n");

}

void isax_topk_query_binary_file(const char *ifilename, int q_num, isax_index *index,
                                 float minimum_distance, int min_checked_leaves, int k,
                                 pqueue_bsf (*search_function)(ts_type *, ts_type *, isax_index *, float, int, int)) {
    fprintf(stderr, ">>> Performing queries in file: %s\n", ifilename);

    FILE *ifile;
    ifile = fopen(ifilename, "rb");
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);
    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }

    int q_loaded = 0;
    ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    ts_type *paa = malloc(sizeof(ts_type) * index->settings->n_segments);

    while (q_loaded < q_num) {
        COUNT_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_INPUT_TIME_END
        // Parse ts and make PAA representation
        paa_from_ts(ts, paa, index->settings);
        COUNT_TOTAL_TIME_START
        COUNT_OUTPUT2_TIME_START
        pqueue_bsf result = search_function(ts, paa, index, minimum_distance, min_checked_leaves, k);
        COUNT_OUTPUT2_TIME_END
        COUNT_TOTAL_TIME_END
        for (int i = 0; i < result.k; i++) {
            printf(" the [%d] query [%d] NN is %f at %ld\n", q_loaded, i, result.knn[i], result.position[i]);
        }
        PRINT_STATS(result.knn[result.k - 1])
        fflush(stdout);

        q_loaded++;
    }
    free(paa);
    free(ts);
    fclose(ifile);
    fprintf(stderr, ">>> Finished querying.\n");

}

void isax_knn_query_binary_file(const char *ifilename, const char *labelfilename, int q_num, isax_index *index,
                                float minimum_distance, int min_checked_leaves, int k, long int classlength,
                                pqueue_bsf (*search_function)(ts_type *, ts_type *, isax_index *, float, int, int)) {
    fprintf(stderr, ">>> Performing queries in file: %s and label in file %s\n", ifilename, labelfilename);

    FILE *ifile, *lfile;
    ifile = fopen(ifilename, "rb");
    lfile = fopen(labelfilename, "rb");
    long int *datalabel;
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    if (lfile == NULL) {
        fprintf(stderr, "File %s not found!\n", labelfilename);
        exit(-1);
    }

    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }
    datalabel = malloc(sizeof(long int) * index->sax_cache_size);
    fread(datalabel, sizeof(long int), index->sax_cache_size, lfile);


    int q_loaded = 0;
    ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    ts_type *paa = malloc(sizeof(ts_type) * index->settings->n_segments);

    while (q_loaded < q_num) {
        COUNT_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_INPUT_TIME_END
        // Parse ts and make PAA representation
        paa_from_ts(ts, paa, index->settings);
        COUNT_TOTAL_TIME_START
        COUNT_OUTPUT2_TIME_START
        pqueue_bsf result = search_function(ts, paa, index, minimum_distance, min_checked_leaves, k);
        COUNT_OUTPUT2_TIME_END
        COUNT_TOTAL_TIME_END
        long int *classcounter = malloc(sizeof(long int) * classlength);
        long int classeposition = 0, classtmp = 0;
        for (int i = 0; i < classlength; i++) {
            classcounter[i] = 0;
        }
        for (int i = 0; i < result.k; i++) {
            classcounter[datalabel[result.position[i]]]++;
            if (classtmp < classcounter[datalabel[result.position[i]]]) {
                classtmp = classcounter[datalabel[result.position[i]]];
                classeposition = datalabel[result.position[i]];

            }
        }
        printf(" the [%d] query's label is %ld \n", q_loaded, classeposition);

        fflush(stdout);

        q_loaded++;
        free(classcounter);
    }
    free(paa);
    free(ts);
    free(datalabel);
    fclose(ifile);
    fclose(lfile);
    fprintf(stderr, ">>> Finished querying.\n");

}


void
isax_knn_query_binary_file_traditional(const char *ifilename, const char *labelfilename, int q_num, isax_index *index,
                                       float minimum_distance, int min_checked_leaves, int k, long int classlength,
                                       pqueue_bsf (*search_function)(ts_type *, ts_type *, isax_index *, node_list *,
                                                                     float, int, int)) {
    fprintf(stderr, ">>> Performing queries in file: %s and label in file %s\n", ifilename, labelfilename);

    FILE *ifile, *lfile;
    ifile = fopen(ifilename, "rb");
    lfile = fopen(labelfilename, "rb");
    long int *datalabel;
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    if (lfile == NULL) {
        fprintf(stderr, "File %s not found!\n", labelfilename);
        exit(-1);
    }

    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }
    datalabel = malloc(sizeof(long int) * index->sax_cache_size);
    fread(datalabel, sizeof(long int), index->sax_cache_size, lfile);


    int q_loaded = 0;
    ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    ts_type *paa = malloc(sizeof(ts_type) * index->settings->n_segments);
    node_list nodelist;
    nodelist.nlist = malloc(sizeof(isax_node *) * pow(2, index->settings->n_segments));
    nodelist.node_amount = 0;
    isax_node *current_root_node = index->first_node;
    while (1) {
        if (current_root_node != NULL) {
            nodelist.nlist[nodelist.node_amount] = current_root_node;
            current_root_node = current_root_node->next;
            nodelist.node_amount++;
        } else {
            break;
        }

    }

    while (q_loaded < q_num) {
        COUNT_INPUT_TIME_START
        fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_INPUT_TIME_END
        PRINT_STATS_HEADER();
        printf("%3d: ", q_loaded);

        // Parse ts and make PAA representation
        paa_from_ts(ts, paa, index->settings);
        COUNT_TOTAL_TIME_START
        pqueue_bsf result = search_function(ts, paa, index, &nodelist, minimum_distance, min_checked_leaves, k);
        COUNT_TOTAL_TIME_END
        long int *classcounter = malloc(sizeof(long int) * classlength);
        long int classeposition = 0, classtmp = 0;
        for (int i = 0; i < classlength; i++) {
            classcounter[i] = 0;
        }
        for (int i = 0; i < result.k; i++) {
            classcounter[datalabel[result.position[i]]]++;
            if (classtmp < classcounter[datalabel[result.position[i]]]) {
                classtmp = classcounter[datalabel[result.position[i]]];
                classeposition = datalabel[result.position[i]];

            }
        }
        printf(" the [%d] query's label is %ld \n", q_loaded, classeposition);

        fflush(stdout);

        q_loaded++;
        free(classcounter);
    }
    free(paa);
    free(ts);
    free(datalabel);
    free(nodelist.nlist);
    fclose(ifile);
    fclose(lfile);
    fprintf(stderr, ">>> Finished querying.\n");

}


void isax_index_binary_file(const char *ifilename, int ts_num, isax_index *index) {
    fprintf(stderr, ">>> Indexing: %s\n", ifilename);
    FILE *ifile;
    COUNT_INPUT_TIME_START
    ifile = fopen(ifilename, "rb");
    COUNT_INPUT_TIME_END
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);
    if (total_records < ts_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }


    int ts_loaded = 0;

    ts_type *ts = malloc(sizeof(ts_type) * index->settings->timeseries_size);
    sax_type *sax = malloc(sizeof(sax_type) * index->settings->n_segments);
    file_position_type *pos = malloc(sizeof(file_position_type));

    index->settings->raw_filename = malloc(256);
    strcpy(index->settings->raw_filename, ifilename);

#ifdef BENCHMARK
    int percentage = (int) (ts_num / (file_position_type) 100);
#endif

    while (ts_loaded < ts_num) {
#ifndef DEBUG
#if VERBOSE_LEVEL == 2
        printf("\r\x1b[32mLoading: \x1b[36m%d\x1b[0m",(ts_loaded + 1));
#endif
#endif
        *pos = ftell(ifile);
        COUNT_INPUT_TIME_START
        int read_number = fread(ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_INPUT_TIME_END

        if (sax_from_ts(ts, sax, index->settings) == SUCCESS) {
#ifdef CLUSTERED
            root_mask_type first_bit_mask = 0x00;
            CREATE_MASK(first_bit_mask, index, sax);
            char* pfilename = malloc(255);
            snprintf(pfilename, 255, "%s.%llu",index->settings->raw_filename,first_bit_mask);
            FILE *pfile = fopen(pfilename, "a+");
            *pos = ftell(pfile);
            fwrite(ts, sizeof(ts_type), index->settings->timeseries_size, pfile);
            fclose(pfile);
            free(pfilename);
#endif
            isax_fbl_index_insert(index, sax, pos);
            ts_loaded++;

        } else {
            fprintf(stderr, "error: cannot insert record in index, since sax representation\
                    failed to be created");
        }
    }
    free(ts);
    free(sax);
    free(pos);
    COUNT_INPUT_TIME_START
    fclose(ifile);
    COUNT_INPUT_TIME_END
    fprintf(stderr, ">>> Finished indexing\n");
}
