#ifdef VALUES
#include <values.h>
#endif
#include <float.h>
#include "../../config.h"
#include "../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>

#include "ads/isax_query_engine.h"
#include "ads/parallel_query_engine.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/pqueue.h"
#include "ads/sax/sax.h"
#include "ads/isax_node_split.h"
#define NTHREADS 4

void isax_query_binary_file_para(const char *ifilename, int q_num, isax_index *index,
        float minimum_distance, int min_checked_leaves,
        query_result (*search_function)(ts_type*, ts_type*, isax_index*, float, int)) {
    fprintf(stderr, ">>> Performing queries in file: %s\n", ifilename);

    FILE * ifile;
    ifile = fopen(ifilename, "rb");

    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);

        exit(-1);
    }

    pthread_t threadid[q_num];
    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz / index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);

    if (total_records < q_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);

        exit(-1);
    }

    int q_loaded = 0;

    paraquery paraqueries[q_num];
    pthread_mutex_t lock_index = PTHREAD_MUTEX_INITIALIZER;
    COUNT_TOTAL_TIME_START
    COUNT_OUTPUT2_TIME_START

    while (q_loaded < q_num) {
        // Parse ts and make PAA representation
        paraqueries[q_loaded].ts = (ts_type*)malloc(sizeof(ts_type) * index->settings->timeseries_size);
        paraqueries[q_loaded].paa = (ts_type*)malloc(sizeof(ts_type) * index->settings->paa_segments);
        paraqueries[q_loaded].index = index;
        paraqueries[q_loaded].minimum_distance = minimum_distance;
        paraqueries[q_loaded].min_checked_leaves = min_checked_leaves;
        paraqueries[q_loaded].lock_index = &lock_index;
        COUNT_INPUT_TIME_START
        fread(paraqueries[q_loaded].ts, sizeof(ts_type), index->settings->timeseries_size, ifile);
        COUNT_INPUT_TIME_END
        paa_from_ts(paraqueries[q_loaded].ts, paraqueries[q_loaded].paa, index->settings->paa_segments,
                index->settings->ts_values_per_paa_segment, index->settings->timeseries_size);
        pthread_create(&(threadid[q_loaded]), NULL, para_queries_worker, (void*)&(paraqueries[q_loaded]));
        fflush(stdout);

#if VERBOSE_LEVEL >= 1
        printf("[%p]: Distance: %lf\n", result.node, result.distance);
#endif

        q_loaded++;
    }

    q_loaded = 0;

    while (q_loaded < q_num) {
        pthread_join(threadid[q_loaded], NULL);
        q_loaded++;
    }

    COUNT_OUTPUT2_TIME_END
    COUNT_TOTAL_TIME_END
    fclose(ifile);
    fprintf(stderr, ">>> Finished querying.\n");
}

void* para_queries_worker(void *transvector) {
    query_result result = exact_search_serial_para(((paraquery*)transvector)->ts, ((paraquery*)transvector)->paa,
            ((paraquery*)transvector)->index, ((paraquery*)transvector)->minimum_distance,
            ((paraquery*)transvector)->min_checked_leaves, ((paraquery*)transvector)->lock_index);
    PRINT_STATS(result.distance);
}

query_result exact_search_serial_para(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance,
        int min_checked_leaves, pthread_mutex_t *lock_index) {
    RESET_BYTES_ACCESSED

    // FOR THREAD USE
    float *MINDISTS = (float*)malloc(sizeof(float) * index->sax_cache_size);

    for (unsigned long j = 0; j < index->sax_cache_size; j++) {
        MINDISTS[j] = FLT_MAX;
    }
    // END

    pthread_mutex_lock(lock_index);
    query_result approximate_result = approximate_search(ts, paa, index);
    query_result bsf_result = approximate_result;

    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }

    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }

    pthread_mutex_unlock(lock_index);

    unsigned long i;
    COUNT_INPUT_TIME_START
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    COUNT_INPUT_TIME_END

    for (i = 0; i < index->sax_cache_size; i++) {
        sax_type *sax = &(index->sax_cache[i * index->settings->paa_segments]);
        MINDISTS[i] = minidist_paa_to_isax_raw(paa, sax,
                index->settings->max_sax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments, MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
    }
    // END

    for (i = 0; i < index->sax_cache_size; i++) {
        sax_type *sax = &index->sax_cache[i * index->settings->paa_segments];

        if (MINDISTS[i] <= approximate_result.distance) {
            COUNT_INPUT_TIME_START
            fseek(raw_file, i * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
            COUNT_INPUT_TIME_END
            float dist = ts_euclidean_distance(ts, ts_buffer, index->settings->timeseries_size, approximate_result.distance);

            if (dist < approximate_result.distance) {
                approximate_result.distance = dist;

#ifdef STORE_ANSWER
                memcpy(index->answer, ts_buffer, index->settings->timeseries_size * sizeof(ts_type));
#endif
            }
        }
    }

    free(ts_buffer);
    fclose(raw_file);
    free(MINDISTS);

    return approximate_result;
}

query_result exact_search_serial_ParIS(ts_type *ts, ts_type *paa, isax_index *index,
        float minimum_distance, int min_checked_leaves) {
    RESET_BYTES_ACCESSED
    pthread_t threadid[maxquerythread];
    query_result approximate_result = approximate_search(ts, paa, index);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    int sum_of_lab = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    float bsf_distance;
    unsigned long j;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }

    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }

    query_result bsf_result = approximate_result;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    unsigned long i;

#ifdef AUTO_TUNE
    float *mindists = malloc(sizeof(float) * index->sax_cache_size);
#endif

    COUNT_CAL_TIME_START
    SET_APPROXIMATE(approximate_result.distance);
    ParIS_LDCW_data *essdata = (ParIS_LDCW_data*)malloc(sizeof(ParIS_LDCW_data) * (maxquerythread));
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;

    for (i = 0; i < (maxquerythread - 1); i++) {
        essdata[i].index = index;
        essdata[i].lock_bsf = &lock_bsf;
        essdata[i].start_number = i * (index->sax_cache_size / maxquerythread);
        essdata[i].stop_number = (i + 1) * (index->sax_cache_size / maxquerythread);
        essdata[i].paa = paa;
        essdata[i].ts = ts;
        essdata[i].bsfdistance = approximate_result.distance;
        essdata[i].sum_of_lab = 0;
    }

    essdata[maxquerythread - 1].index = index;
    essdata[maxquerythread - 1].lock_bsf = &lock_bsf;
    essdata[maxquerythread - 1].start_number = (maxquerythread - 1) * (index->sax_cache_size / maxquerythread);
    essdata[maxquerythread - 1].stop_number = index->sax_cache_size;
    essdata[maxquerythread - 1].paa = paa;
    essdata[maxquerythread - 1].ts = ts;
    essdata[maxquerythread - 1].bsfdistance = approximate_result.distance;
    essdata[maxquerythread - 1].sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, mindistance_worker, (void*)&(essdata[i]));
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
        sum_of_lab += essdata[i].sum_of_lab;
    }

    unsigned long* label_number = (unsigned long*)malloc(sizeof(unsigned long) * (sum_of_lab));
    float* minidisvector = (float*)malloc(sizeof(float) * (sum_of_lab));

    sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        memcpy(&(label_number[sum_of_lab]), essdata[i].label_number, sizeof(unsigned long) * essdata[i].sum_of_lab);
        memcpy(&(minidisvector[sum_of_lab]), essdata[i].minidisvector, sizeof(float) * essdata[i].sum_of_lab);
        free(essdata[i].label_number);
        free(essdata[i].minidisvector);
        sum_of_lab += essdata[i].sum_of_lab;
    }

    pthread_t readthread[maxquerythread * maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter = 0;
    float bsfdistance = (approximate_result.distance);
    readpointer.ts = ts;
    readpointer.index = index;
    readpointer.counter = &readcounter;
    readpointer.bsf = approximate_result.distance;
    readpointer.load_point = label_number;
    readpointer.lock_bsf = &lock_bsf;
    readpointer.bsf2 = &bsfdistance;
    readpointer.minidisvector = minidisvector;
    readpointer.sum_of_lab = sum_of_lab;
    COUNT_CAL_TIME_END
    COUNT_INPUT_TIME_START

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_create(&(readthread[i]), NULL, read_worker, (void*)&(readpointer));
    }

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_join(readthread[i], NULL);
    }

    COUNT_INPUT_TIME_END
    approximate_result.distance = bsfdistance;
    free(essdata);
    free(minidisvector);
    free(label_number);
    fclose(raw_file);

    return approximate_result;
}

query_result exact_search_serial_ParISnonsort(ts_type *ts, ts_type *paa, isax_index *index,
        float minimum_distance, int min_checked_leaves) {
    RESET_BYTES_ACCESSED
    pthread_t threadid[maxquerythread];
    query_result approximate_result = approximate_search(ts, paa, index);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    unsigned long sum_of_lab = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    float bsf_distance;
    unsigned long j;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }

    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }

    query_result bsf_result = approximate_result;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    unsigned long i;

#ifdef AUTO_TUNE
    float *mindists = (float*)malloc(sizeof(float) * index->sax_cache_size);
#endif

    COUNT_CAL_TIME_START
    SET_APPROXIMATE(approximate_result.distance);
    ParIS_LDCW_data *essdata = (ParIS_LDCW_data*)malloc(sizeof(ParIS_LDCW_data)*(maxquerythread));
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;

    for (i = 0; i < (maxquerythread - 1); i++) {
        essdata[i].index = index;
        essdata[i].lock_bsf = &lock_bsf;
        essdata[i].start_number = i * (index->sax_cache_size / maxquerythread);
        essdata[i].stop_number = (i + 1) * (index->sax_cache_size / maxquerythread);
        essdata[i].paa = paa;
        essdata[i].ts = ts;
        essdata[i].bsfdistance = approximate_result.distance;
        essdata[i].sum_of_lab = 0;
    }

    essdata[maxquerythread - 1].index = index;
    essdata[maxquerythread - 1].lock_bsf = &lock_bsf;
    essdata[maxquerythread - 1].start_number = (maxquerythread - 1) * (index->sax_cache_size / maxquerythread);
    essdata[maxquerythread - 1].stop_number = index->sax_cache_size;
    essdata[maxquerythread - 1].paa = paa;
    essdata[maxquerythread - 1].ts = ts;
    essdata[maxquerythread - 1].bsfdistance = approximate_result.distance;
    essdata[maxquerythread - 1].sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, mindistance_worker, (void*)&(essdata[i]));
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
        sum_of_lab += essdata[i].sum_of_lab;
    }

    unsigned long* label_number = (unsigned long*)malloc(sizeof(unsigned long) * (sum_of_lab));
    float* minidisvector = (float*)malloc(sizeof(float) * (sum_of_lab));
    sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        essdata[i].currentpositioncounter = &sum_of_lab;
        essdata[i].label_number = label_number;
        essdata[i].minidisvector = minidisvector;
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, mindistanceinsert_worker, (void*)&(essdata[i]));
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    pthread_t readthread[maxquerythread * maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter = 0;
    float bsfdistance = (approximate_result.distance);

    readpointer.ts = ts;
    readpointer.index = index;
    readpointer.counter = &readcounter;
    readpointer.bsf = approximate_result.distance;
    readpointer.load_point = label_number;
    readpointer.lock_bsf = &lock_bsf;
    readpointer.bsf2 = &bsfdistance;
    readpointer.minidisvector = minidisvector;
    readpointer.sum_of_lab = sum_of_lab;
    COUNT_CAL_TIME_END
    COUNT_INPUT_TIME_START

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_create(&(readthread[i]), NULL, read_worker, (void*)&(readpointer));
    }

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_join(readthread[i], NULL);
    }

    COUNT_INPUT_TIME_END
    approximate_result.distance = bsfdistance;
    free(essdata);
    free(minidisvector);
    free(label_number);
    fclose(raw_file);

    return approximate_result;
}

pqueue_bsf exact_topk_serial_ParIS(ts_type *ts, ts_type *paa, isax_index *index,
        float minimum_distance, int min_checked_leaves, int k) {
    RESET_BYTES_ACCESSED
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    pthread_t threadid[maxquerythread];
    pqueue_bsf *pq_bsf = pqueue_bsf_init(k);
    approximate_topk(ts, paa, index, pq_bsf);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    int sum_of_lab = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    float bsf_distance;
    unsigned long j;

    if (pq_bsf->knn[k - 1] == 0) {
        return *pq_bsf;
    }

    if (pq_bsf->knn[k - 1] == FLT_MAX  || min_checked_leaves > 1) {
        refine_topk_answer(ts, paa, index, pq_bsf, minimum_distance, min_checked_leaves);
    }

    unsigned long i;

#ifdef AUTO_TUNE
    float *mindists = (float*)malloc(sizeof(float) * index->sax_cache_size);
#endif

    SET_APPROXIMATE(pq_bsf->knn[k - 1]);
    ParIS_LDCW_data *essdata = (ParIS_LDCW_data*)malloc(sizeof(ParIS_LDCW_data) * (maxquerythread));
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;

    for (i = 0; i < (maxquerythread - 1); i++) {
        essdata[i].index = index;
        essdata[i].lock_bsf = &lock_bsf;
        essdata[i].start_number = i * (index->sax_cache_size / maxquerythread);
        essdata[i].stop_number = (i + 1) * (index->sax_cache_size / maxquerythread);
        essdata[i].paa = paa;
        essdata[i].ts = ts;
        essdata[i].bsfdistance = pq_bsf->knn[k - 1];
        essdata[i].sum_of_lab = 0;
    }

    essdata[maxquerythread - 1].index = index;
    essdata[maxquerythread - 1].lock_bsf = &lock_bsf;
    essdata[maxquerythread - 1].start_number = (maxquerythread - 1) * (index->sax_cache_size / maxquerythread);
    essdata[maxquerythread - 1].stop_number = index->sax_cache_size;
    essdata[maxquerythread - 1].paa = paa;
    essdata[maxquerythread - 1].ts = ts;
    essdata[maxquerythread - 1].bsfdistance = pq_bsf->knn[k - 1];
    essdata[maxquerythread - 1].sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, mindistance_worker, (void*)&(essdata[i]));
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
        sum_of_lab += essdata[i].sum_of_lab;
    }

    unsigned long* label_number = (unsigned long*)malloc(sizeof(unsigned long) * (sum_of_lab));
    float* minidisvector = (float*)malloc(sizeof(float) * (sum_of_lab));
    sum_of_lab = 0;

    for (i = 0; i < maxquerythread; i++) {
        memcpy(&(label_number[sum_of_lab]), essdata[i].label_number, sizeof(unsigned long) * essdata[i].sum_of_lab);
        memcpy(&(minidisvector[sum_of_lab]), essdata[i].minidisvector, sizeof(float) * essdata[i].sum_of_lab);
        free(essdata[i].label_number);
        free(essdata[i].minidisvector);
        sum_of_lab += essdata[i].sum_of_lab;
    }

    pthread_t readthread[maxquerythread * maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter = 0;
    readpointer.ts = ts;
    readpointer.index = index;
    readpointer.counter = &readcounter;
    readpointer.load_point = label_number;
    readpointer.lock_bsf = &lock_bsf;
    readpointer.minidisvector = minidisvector;
    readpointer.sum_of_lab = sum_of_lab;
    readpointer.pq_bsf = pq_bsf;
    unsigned long read_time_conter[maxquerythread * maxreadthread];
    unsigned long read_time_all = 0;

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_create(&(readthread[i]), NULL, topk_read_worker, (void*)&(readpointer));
    }

    for (i = 0; i < maxquerythread * maxreadthread; i++) {
        pthread_join(readthread[i], NULL);
    }

    free(essdata);
    free(minidisvector);
    free(label_number);
    fclose(raw_file);

    return *pq_bsf;
}

void* mindistanceinsert_worker(void *essdata) {
    unsigned long* localposition = ((ParIS_LDCW_data*)essdata)->label_number;
    float* localmindist = ((ParIS_LDCW_data*)essdata)->minidisvector;
    isax_index *index = ((ParIS_LDCW_data*)essdata)->index;
    unsigned long start_number = ((ParIS_LDCW_data*)essdata)->start_number;
    unsigned long stop_number = ((ParIS_LDCW_data*)essdata)->stop_number;
    unsigned long t;
    float bsfdistance;
    float mindist;
    ts_type *paa = ((ParIS_LDCW_data*)essdata)->paa;
    ts_type *ts = ((ParIS_LDCW_data*)essdata)->ts;

    for (unsigned long i = start_number; i < stop_number; i++) {
        sax_type *sax = &index->sax_cache[i * index->settings->paa_segments];
#if defined(__x86_64__)
        mindist = minidist_paa_to_isax_rawa_SIMD(paa, sax, index->settings->max_sax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments, MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
#else
        mindist = minidist_paa_to_isax_raw(paa, sax, index->settings->max_sax_cardinalities,
                 index->settings->sax_bit_cardinality,
                 index->settings->sax_alphabet_cardinality,
                 index->settings->paa_segments, MINVAL, MAXVAL,
                 index->settings->mindist_sqrt);
#endif
        if (mindist <= ((ParIS_LDCW_data*)essdata)->bsfdistance) {
            t = __sync_fetch_and_add(((ParIS_LDCW_data*)essdata)->currentpositioncounter, 1);
            memcpy(&(((ParIS_LDCW_data*)essdata)->label_number[t]), &i, sizeof(unsigned int));
            memcpy(&(((ParIS_LDCW_data*)essdata)->minidisvector[t]), &mindist, sizeof(float));
        }
    }
}

void* mindistance_worker(void *essdata) {
    float bsfdistance;
    float mindist;
    isax_index *index = ((ParIS_LDCW_data*)essdata)->index;
    unsigned long start_number = ((ParIS_LDCW_data*)essdata)->start_number;
    unsigned long stop_number = ((ParIS_LDCW_data*)essdata)->stop_number;
    ts_type *paa = ((ParIS_LDCW_data*)essdata)->paa;
    ts_type *ts = ((ParIS_LDCW_data*)essdata)->ts;
    ((ParIS_LDCW_data*)essdata)->label_number = (unsigned long*)malloc(sizeof(unsigned long) * 10000);
    ((ParIS_LDCW_data*)essdata)->minidisvector = (float*)malloc(sizeof(float) * 10000);
    unsigned long max_number = 10000;

    for (unsigned long i = start_number; i < stop_number; i++) {
        sax_type *sax = &index->sax_cache[i * index->settings->paa_segments];
#if defined(__x86_64__)
        mindist = minidist_paa_to_isax_rawa_SIMD(paa, sax, index->settings->max_sax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments, MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
#else
        mindist = minidist_paa_to_isax_raw(paa, sax, index->settings->max_sax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments, MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
#endif

        if (mindist <= ((ParIS_LDCW_data*)essdata)->bsfdistance) {
            if (((ParIS_LDCW_data*)essdata)->sum_of_lab >= max_number) {
                max_number = (max_number + 10000);
                unsigned long* change_lab = ((ParIS_LDCW_data*)essdata)->label_number;
                float* change_minivec = ((ParIS_LDCW_data*)essdata)->minidisvector;
                ((ParIS_LDCW_data*)essdata)->label_number = (unsigned long*)malloc(sizeof(unsigned long) * (max_number + 10000));
                ((ParIS_LDCW_data*)essdata)->minidisvector = (float*)malloc(sizeof(float) * (max_number + 10000));
                memcpy(((ParIS_LDCW_data*)essdata)->label_number, change_lab, sizeof(unsigned long) * max_number);
                memcpy(((ParIS_LDCW_data*)essdata)->minidisvector, change_minivec, sizeof(float) * max_number);
                free(change_lab);
                free(change_minivec);
            }

            ((ParIS_LDCW_data*)essdata)->label_number[((ParIS_LDCW_data*)essdata)->sum_of_lab] = i;
            ((ParIS_LDCW_data*)essdata)->minidisvector[((ParIS_LDCW_data*)essdata)->sum_of_lab] = mindist;
            ((ParIS_LDCW_data*)essdata)->sum_of_lab++;
        }
    }
}

void* read_worker(void *read_pointer) {
    isax_index *index = ((ParIS_read_worker_data*)read_pointer)->index;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    ts_type *ts =((ParIS_read_worker_data*)read_pointer)->ts;
    unsigned long t = 0;
    unsigned long p;
    unsigned long sum_of_lab = ((ParIS_read_worker_data*)read_pointer)->sum_of_lab;
    float *minidisvector = ((ParIS_read_worker_data*)read_pointer)->minidisvector;
    float bsf;
    float dist;

    while (1) {
        pthread_rwlock_rdlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf);
        bsf = *(((ParIS_read_worker_data*)read_pointer)->bsf2);
        pthread_rwlock_unlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf);
        t = __sync_fetch_and_add(((ParIS_read_worker_data*)read_pointer)->counter, 1);

        if (t >= sum_of_lab) {
            break;
        }

        p = ((ParIS_read_worker_data*)read_pointer)->load_point[t];

        if (minidisvector[t] < bsf) {
            fseek(raw_file, p * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
#if defined(__x86_64__)
            dist = ts_euclidean_distance_SIMD(ts, ts_buffer, index->settings->timeseries_size, bsf);
#else
            dist = ts_euclidean_distance(ts, ts_buffer, index->settings->timeseries_size, bsf);
#endif

            if (dist < bsf) {
                pthread_rwlock_wrlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf);

                if (dist < *(((ParIS_read_worker_data*)read_pointer)->bsf2)) {
                    *(((ParIS_read_worker_data*)read_pointer)->bsf2) = dist;
                }

                pthread_rwlock_unlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf);
            }
        }
    }

    free(ts_buffer);
    fclose(raw_file);
}

void* topk_read_worker(void *read_pointer) {
    isax_index *index = ((ParIS_read_worker_data*)read_pointer)->index;
    pqueue_bsf *pq_bsf = ((ParIS_read_worker_data*)read_pointer)->pq_bsf;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    ts_type *ts = ((ParIS_read_worker_data*)read_pointer)->ts;
    unsigned long t = 0;
    unsigned long p;
    unsigned long sum_of_lab = ((ParIS_read_worker_data*)read_pointer)->sum_of_lab;
    float *minidisvector = ((ParIS_read_worker_data*)read_pointer)->minidisvector;
    float bsf;
    float dist;

    while (1) {
        bsf = pq_bsf->knn[pq_bsf->k - 1];
        t = __sync_fetch_and_add(((ParIS_read_worker_data*)read_pointer)->counter, 1);

        if (t >= sum_of_lab) {
            break;
        }

        p = ((ParIS_read_worker_data*)read_pointer)->load_point[t];

        if (minidisvector[t] < bsf) {
            fseek(raw_file, p * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
#if defined(__x86_64__)
            dist = ts_euclidean_distance_SIMD(ts, ts_buffer, index->settings->timeseries_size, bsf);
#else
            dist = ts_euclidean_distance(ts, ts_buffer, index->settings->timeseries_size, bsf);
#endif

            if (dist <= bsf) {
                pthread_rwlock_wrlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf);
                pqueue_bsf_insert(pq_bsf, dist, p, NULL);
                pthread_rwlock_unlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf);
            }
        }
    }

    free(ts_buffer);
    fclose(raw_file);
}

query_result exact_search_serial_ParIS_nb(ts_type *ts, ts_type *paa, isax_index *index,
        float minimum_distance, int min_checked_leaves)  {
    RESET_BYTES_ACCESSED
    pthread_t threadid[maxquerythread];
    query_result approximate_result = approximate_search(ts, paa, index);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    query_result bsf_result = approximate_result;
    int sum_of_lab = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    float bsf_distance;
    unsigned long j;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }

    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }

    unsigned long i;

#ifdef AUTO_TUNE
    float *mindists = (float*)malloc(sizeof(float) * index->sax_cache_size);
#endif

    SET_APPROXIMATE(approximate_result.distance);
    ParIS_LDCW_data *essdata = (ParIS_LDCW_data*)malloc(sizeof(ParIS_LDCW_data) * (maxquerythread));
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    float bsfdistance = approximate_result.distance;
    unsigned long read_time_conter = 0;

    for (i = 0; i < (maxquerythread - 1); i++) {
        essdata[i].index = index;
        essdata[i].lock_bsf = &lock_bsf;
        essdata[i].start_number = i * (index->sax_cache_size / maxquerythread);
        essdata[i].stop_number = (i + 1) * (index->sax_cache_size / maxquerythread);
        essdata[i].paa = paa;
        essdata[i].ts = ts;
        essdata[i].bsfdistance = approximate_result.distance;
        essdata[i].sum_of_lab = 0;
        essdata[i].minidisvector = &bsfdistance;
        essdata[i].read_time_conter = 0;
    }

    essdata[maxquerythread - 1].index = index;
    essdata[maxquerythread - 1].lock_bsf = &lock_bsf;
    essdata[maxquerythread - 1].start_number = (maxquerythread - 1) * (index->sax_cache_size / maxquerythread);
    essdata[maxquerythread - 1].stop_number = index->sax_cache_size;
    essdata[maxquerythread - 1].paa = paa;
    essdata[maxquerythread - 1].ts = ts;
    essdata[maxquerythread - 1].bsfdistance = approximate_result.distance;
    essdata[maxquerythread - 1].sum_of_lab = 0;
    essdata[maxquerythread - 1].minidisvector = &bsfdistance;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, ParIS_nb_worker, (void*)&(essdata[i]));
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);

        if (essdata[i].bsfdistance < approximate_result.distance) {
            approximate_result.distance = essdata[i].bsfdistance;
        }
    }

    free(essdata);

    return approximate_result;
}

void* ParIS_nb_worker(void *essdata) {
    unsigned long i;
    float bsfdistance;
    float mindist;
    isax_index *index = ((ParIS_LDCW_data*)essdata)->index;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = (ts_type*)malloc(index->settings->ts_byte_size);
    unsigned long start_number = ((ParIS_LDCW_data*)essdata)->start_number;
    unsigned long stop_number = ((ParIS_LDCW_data*)essdata)->stop_number;
    ts_type *paa = ((ParIS_LDCW_data*)essdata)->paa;
    ts_type *ts = ((ParIS_LDCW_data*)essdata)->ts;

    for (i = start_number; i < stop_number; i++) {
        sax_type *sax = &index->sax_cache[i * index->settings->paa_segments];
#if defined(__x86_64__)
        mindist = minidist_paa_to_isax_raw_SIMD(paa, sax, index->settings->max_sax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments, MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
#else
        mindist = minidist_paa_to_isax_raw(paa, sax, index->settings->max_sax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments, MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
#endif

        if (mindist <= ((ParIS_LDCW_data*)essdata)->bsfdistance) {
            fseek(raw_file, i * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
#if defined(__x86_64__)
            float dist = ts_euclidean_distance_SIMD(ts, ts_buffer,
                    index->settings->timeseries_size, ((ParIS_LDCW_data*)essdata)->bsfdistance);
#else
            float dist = ts_euclidean_distance(ts, ts_buffer,
                    index->settings->timeseries_size, ((ParIS_LDCW_data*)essdata)->bsfdistance);
#endif

            if (dist < (((ParIS_LDCW_data*)essdata)->bsfdistance)) {
                (((ParIS_LDCW_data*)essdata)->bsfdistance) = dist;
            }
        }
    }

    fclose(raw_file);
    free(ts_buffer);
}

query_result refin_answer_m(ts_type *ts, ts_type *paa, isax_index *index,
        query_result *approximate_bsf_result, float minimum_distance, int limit) {
    int maxid = 8;
    pthread_t threadid[maxid];
    refind_answer_fonction_data rfdata;
    pthread_mutex_t lock_queue = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t lock_current_root_node = PTHREAD_MUTEX_INITIALIZER;
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    pthread_barrier_t lock_barrier;
    pthread_barrier_init(&lock_barrier, NULL, 4);
    rfdata.paa = paa;
    rfdata.ts = ts;
    rfdata.lock_queue = &lock_queue;
    rfdata.lock_current_root_node = &lock_current_root_node;
    rfdata.lock_bsf = &lock_bsf;
    rfdata.index = index;
    rfdata.minimum_distance = minimum_distance;
    rfdata.limit = limit / 4;
    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size,
            cmp_pri, get_pri, set_pri, get_pos, set_pos);
    rfdata.pq = pq;

    // Insert all root nodes in heap.
    rfdata.current_root_node = index->first_node;
    rfdata.bsf_result = approximate_bsf_result;
    query_result bsf_result;
    query_result *n;

    rfdata.lock_barrier = &lock_barrier;

    for (int i = 0; i < maxid; i++) {
        pthread_create(&(threadid[i]), NULL, refind_answer_fonction, (void*)&(rfdata));
    }

    for (int i = 0; i < maxid; i++) {
        pthread_join(threadid[i], NULL);
    }

    // Free the nodes that where not popped.
    while ((n = (query_result*)pqueue_pop(pq))) {
        free(n);
    }

    // Free the priority queue.
    pthread_barrier_destroy(&lock_barrier);
    pqueue_free(pq);
    bsf_result = *(rfdata.bsf_result);

    return *(rfdata.bsf_result);
}

void* refind_answer_fonction(void *rfdata) {
    isax_node *current_root_node;
    query_result *n;
    isax_index *index = ((refind_answer_fonction_data*)rfdata)->index;
    ts_type *paa = ((refind_answer_fonction_data*)rfdata)->paa;
    ts_type *ts = ((refind_answer_fonction_data*)rfdata)->ts;
    pqueue_t *pq = ((refind_answer_fonction_data*)rfdata)->pq;
    float minimum_distance = ((refind_answer_fonction_data*)rfdata)->minimum_distance;
    float bsfdisntance;
    int limit = ((refind_answer_fonction_data*)rfdata)->limit;
    int checks = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    query_result *bsf_result = (((refind_answer_fonction_data*)rfdata)->bsf_result);

    while (1) {
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
        current_root_node = ((refind_answer_fonction_data*)rfdata)->current_root_node;

        if (current_root_node != NULL) {
            ((refind_answer_fonction_data*)rfdata)->current_root_node =
                    ((refind_answer_fonction_data*)rfdata)->current_root_node->next;
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
        } else {
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);

            break;
        }

        query_result * mindist_result = (query_result*)malloc(sizeof(query_result));
        mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values,
                current_root_node->isax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments,
                MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
        mindist_result->node = current_root_node;
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);

        pqueue_insert(pq, mindist_result);
        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
    }

    pthread_barrier_wait(((refind_answer_fonction_data*)rfdata)->lock_barrier);

    while (1) {
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
        n = (query_result*)pqueue_pop(pq);
        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

        pthread_rwlock_rdlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
        bsfdisntance = bsf_result->distance;
        pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

        // The best node has a worse mindist, so search is finished!
        if (n->distance >= bsfdisntance || n->distance > minimum_distance) {
            pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
            pqueue_insert(pq, n);
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

            break;
        } else {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                // *** ADAPTIVE SPLITTING ***
                if (!n->node->has_full_data_file && (n->node->leaf_size > index->settings->min_leaf_size)) {
                    // Split and push again in the queue
                    split_node(index, n->node);
                    pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    pqueue_insert(pq, n);
                    pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

                    continue;
                }

                // *** EXTRA BOUNDING ***
                if (tight_bound) {
                    float mindistance = calculate_minimum_distance(index, n->node, ts, paa);

                    if (mindistance >= bsfdisntance) {
                        free(n);

                        continue;
                    }
                }

                // *** REAL DISTANCE ***
                checks++;
                float distance = calculate_node_distance(index, n->node, ts, bsfdisntance);

                if (distance < bsfdisntance) {
                    pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

                    if (distance < bsf_result->distance) {
                        bsf_result->distance = distance;
                        bsf_result->node = n->node;
                    }

                    pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                }

                if (checks > limit) {
                    pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    pqueue_insert(pq, n);
                    pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

                    break;
                }
            } else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if (n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check) {
                        float distance = calculate_node_distance(index, n->node->left_child, ts, bsfdisntance);

                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

                            if (distance < bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->left_child;
                            }

                            pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                        }
                    } else {
                        query_result *mindist_result = (query_result*)malloc(sizeof(query_result));
                        mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values,
                                n->node->left_child->isax_cardinalities,
                                index->settings->sax_bit_cardinality,
                                index->settings->sax_alphabet_cardinality,
                                index->settings->paa_segments,
                                MINVAL, MAXVAL,
                                index->settings->mindist_sqrt);
                        mindist_result->node = n->node->left_child;
                        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                        pqueue_insert(pq, mindist_result);
                        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    }
                }

                if (n->node->right_child->isax_cardinalities != NULL) {
                    if (n->node->right_child->is_leaf &&
                            !n->node->left_child->has_partial_data_file && aggressive_check) {
                        float distance = calculate_node_distance(index, n->node->right_child, ts, bsfdisntance);

                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

                            if (distance <bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->right_child;
                            }

                            pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                        }
                    } else {
                        query_result *mindist_result = (query_result*)malloc(sizeof(query_result));
                        mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values,
                                n->node->right_child->isax_cardinalities,
                                index->settings->sax_bit_cardinality,
                                index->settings->sax_alphabet_cardinality,
                                index->settings->paa_segments,
                                MINVAL, MAXVAL,
                                index->settings->mindist_sqrt);
                        mindist_result->node = n->node->right_child;
                        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                        pqueue_insert(pq, mindist_result);
                        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    }
                }
            }

            // Free the node currently popped.
            free(n);
        }
    }
}

query_result exact_search_m(ts_type *ts, ts_type *paa, isax_index *index,
        float minimum_distance, int min_checked_leaves) {
    query_result approximate_result = approximate_search_SIMD(ts, paa, index);
    query_result bsf_result = approximate_result;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    int i;

    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }

    if (approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }

    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size,
            cmp_pri, get_pri, set_pri, get_pos, set_pos);
    query_result *do_not_remove = &approximate_result;
    SET_APPROXIMATE(approximate_result.distance);

    RESET_BYTES_ACCESSED

    if (approximate_result.node != NULL) {
        // Insert approximate result in heap.
        pqueue_insert(pq, &approximate_result);
        // GOOD: if(approximate_result.node->filename != NULL)
        // GOOD: printf("POPS: %.5lf\t", approximate_result.distance);
    }

    // Insert all root nodes in heap.
    isax_node *current_root_node = index->first_node;
    pthread_t threadid[maxquerythread];

    refind_answer_fonction_data rfdata;
    pthread_mutex_t lock_queue = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t lock_current_root_node = PTHREAD_MUTEX_INITIALIZER;
    pthread_rwlock_t lock_bsf = PTHREAD_RWLOCK_INITIALIZER;
    pthread_barrier_t lock_barrier;
    pthread_barrier_init(&lock_barrier, NULL, maxquerythread);
    rfdata.paa = paa;
    rfdata.ts = ts;
    rfdata.lock_queue = &lock_queue;
    rfdata.lock_current_root_node = &lock_current_root_node;
    rfdata.lock_bsf = &lock_bsf;
    rfdata.index = index;
    rfdata.minimum_distance = minimum_distance;
    rfdata.pq = pq;

    // Insert all root nodes in heap.
    rfdata.current_root_node = index->first_node;
    rfdata.bsf_result = &bsf_result;

    query_result * n;

    rfdata.lock_barrier = &lock_barrier;

    for (i = 0; i < maxquerythread; i++) {
        pthread_create(&(threadid[i]), NULL, exact_search_fonction, (void*)&(rfdata));
    }

    for (i = 0; i < maxquerythread; i++) {
        pthread_join(threadid[i], NULL);
    }

    // Free the nodes that where not popped.
    while ((n = (query_result*)pqueue_pop(pq))) {
        if (n != do_not_remove) {
            free(n);
        }
    }

    // Free the priority queue.
    pthread_barrier_destroy(&lock_barrier);
    pqueue_free(pq);

    bsf_result = *(rfdata.bsf_result);

    return *(rfdata.bsf_result);
}

void* exact_search_fonction(void *rfdata) {
    isax_node *current_root_node;
    query_result *n;
    isax_index *index = ((refind_answer_fonction_data*)rfdata)->index;
    ts_type *paa = ((refind_answer_fonction_data*)rfdata)->paa;
    ts_type *ts = ((refind_answer_fonction_data*)rfdata)->ts;
    pqueue_t *pq = ((refind_answer_fonction_data*)rfdata)->pq;
    query_result *do_not_remove = ((refind_answer_fonction_data*)rfdata)->bsf_result;
    float minimum_distance = ((refind_answer_fonction_data*)rfdata)->minimum_distance;
    float bsfdisntance;
    int limit = ((refind_answer_fonction_data*)rfdata)->limit;
    int checks = 0;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    query_result *bsf_result = (((refind_answer_fonction_data*)rfdata)->bsf_result);

    while (1) {
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
        current_root_node = ((refind_answer_fonction_data*)rfdata)->current_root_node;

        if (current_root_node != NULL) {
            ((refind_answer_fonction_data*)rfdata)->current_root_node =
                    ((refind_answer_fonction_data*)rfdata)->current_root_node->next;
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
        } else {
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);

            break;
        }

        query_result * mindist_result = (query_result*)malloc(sizeof(query_result));
        mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values,
                current_root_node->isax_cardinalities,
                index->settings->sax_bit_cardinality,
                index->settings->sax_alphabet_cardinality,
                index->settings->paa_segments,
                MINVAL, MAXVAL,
                index->settings->mindist_sqrt);
        mindist_result->node = current_root_node;
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);

        pqueue_insert(pq, mindist_result);

        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
    }

    pthread_barrier_wait(((refind_answer_fonction_data*)rfdata)->lock_barrier);

    while (1) {
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
        n = (query_result*)pqueue_pop(pq);
        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

        if (n == NULL) {
            break;
        }

        pthread_rwlock_rdlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
        bsfdisntance = bsf_result->distance;
        pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

        // The best node has a worse mindist, so search is finished!
        if (n->distance >= bsfdisntance || n->distance > minimum_distance) {
            pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
            pqueue_insert(pq, n);
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

            break;
        } else {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                // *** ADAPTIVE SPLITTING ***
                if (!n->node->has_full_data_file && (n->node->leaf_size > index->settings->min_leaf_size)) {
                    // Split and push again in the queue
                    split_node(index, n->node);
                    pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    pqueue_insert(pq, n);
                    pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);

                    continue;
                }

                // *** EXTRA BOUNDING ***
                if (tight_bound) {
                    float mindistance = calculate_minimum_distance_SIMD(index, n->node, ts, paa);

                    if (mindistance >= bsfdisntance) {
                        if (n != do_not_remove) {
                            free(n);
                        }

                        continue;
                    }
                }

                // *** REAL DISTANCE ***
                checks++;
                float distance = calculate_node_distance_SIMD(index, n->node, ts, bsfdisntance);

                if (distance < bsfdisntance) {
                    pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

                    if (distance < bsf_result->distance) {
                        bsf_result->distance = distance;
                        bsf_result->node = n->node;
                    }

                    pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                }
            } else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if (n->node->left_child->is_leaf &&
                            !n->node->left_child->has_partial_data_file && aggressive_check) {
                        float distance = calculate_node_distance_SIMD(index, n->node->left_child, ts, bsfdisntance);

                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

                            if (distance <bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->left_child;
                            }

                            pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                        }
                    } else {
                        query_result * mindist_result = (query_result*)malloc(sizeof(query_result));
                        mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values,
                                n->node->left_child->isax_cardinalities,
                                index->settings->sax_bit_cardinality,
                                index->settings->sax_alphabet_cardinality,
                                index->settings->paa_segments,
                                MINVAL, MAXVAL,
                                index->settings->mindist_sqrt);
                        mindist_result->node = n->node->left_child;
                        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                        pqueue_insert(pq, mindist_result);
                        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    }
                }

                if (n->node->right_child->isax_cardinalities != NULL) {
                    if (n->node->right_child->is_leaf &&
                            !n->node->left_child->has_partial_data_file && aggressive_check) {
                        float distance = calculate_node_distance_SIMD(index, n->node->right_child, ts, bsfdisntance);

                        if (distance < bsfdisntance) {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);

                            if (distance <bsf_result->distance) {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->right_child;
                            }

                            pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                        }
                    } else {
                        query_result * mindist_result = (query_result*)malloc(sizeof(query_result));
                        mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values,
                                n->node->right_child->isax_cardinalities,
                                index->settings->sax_bit_cardinality,
                                index->settings->sax_alphabet_cardinality,
                                index->settings->paa_segments,
                                MINVAL, MAXVAL,
                                index->settings->mindist_sqrt);
                        mindist_result->node = n->node->right_child;
                        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                        pqueue_insert(pq, mindist_result);
                        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    }
                }
            }

            // Free the node currently popped.
            if (n != do_not_remove) {
                free(n);
            }
        }
    }
}
