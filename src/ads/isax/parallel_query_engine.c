#ifdef VALUES
#include <values.h>
#endif
#include <float.h>
#include "config.h"
#include "../../../globals.h"
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
#include "ads/pthread_barrier.h"

static void *mindistance_worker(void *essdata);
static void *read_worker(void *read_pointer);
static void *topk_read_worker(void *read_pointer);
static void *ParIS_nb_worker(void *essdata);
static void *exact_search_fonction(void *rfdata);

query_result exact_search_serial_ParIS(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves) 
{

    RESET_BYTES_ACCESSED


    pthread_t threadid[maxquerythread];
    query_result approximate_result = approximate_search(ts, paa, index);
    int sum_of_lab=0;
    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    
    if(approximate_result.distance == FLT_MAX || min_checked_leaves > 1)
    {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }
    int i;

    COUNT_CAL_TIME_START
    SET_APPROXIMATE(approximate_result.distance);
    ParIS_LDCW_data *essdata=malloc(sizeof(ParIS_LDCW_data)*(maxquerythread));
    pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;
    for (i = 0; i < (maxquerythread-1); i++)
    {
        essdata[i].index=index;
        essdata[i].lock_bsf=&lock_bsf;
        essdata[i].start_number=i*(index->sax_cache_size/maxquerythread);
        essdata[i].stop_number=(i+1)*(index->sax_cache_size/maxquerythread);
        essdata[i].paa=paa;
        essdata[i].ts=ts;
        essdata[i].bsfdistance=approximate_result.distance;
        essdata[i].sum_of_lab=0;
    }
    essdata[maxquerythread-1].index=index;
    essdata[maxquerythread-1].lock_bsf=&lock_bsf;
    essdata[maxquerythread-1].start_number=(maxquerythread-1)*(index->sax_cache_size/maxquerythread);
    essdata[maxquerythread-1].stop_number=index->sax_cache_size;
    essdata[maxquerythread-1].paa=paa;
    essdata[maxquerythread-1].ts=ts;
    essdata[maxquerythread-1].bsfdistance=approximate_result.distance;
    essdata[maxquerythread-1].sum_of_lab=0;

    for(i=0; i<maxquerythread; i++) 
    {
        pthread_create(&(threadid[i]),NULL,mindistance_worker,(void*)&(essdata[i]));
    }
    for (i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
        sum_of_lab+=essdata[i].sum_of_lab;
    }

    unsigned long* label_number=malloc(sizeof(unsigned long)*(sum_of_lab));
    float* minidisvector=malloc(sizeof(float)*(sum_of_lab));
   
    sum_of_lab=0;

    for (i = 0; i < maxquerythread; i++)
    {
        memcpy(&(label_number[sum_of_lab]),essdata[i].label_number,sizeof(unsigned long)*essdata[i].sum_of_lab);
        memcpy(&(minidisvector[sum_of_lab]),essdata[i].minidisvector,sizeof(float)*essdata[i].sum_of_lab);
        free(essdata[i].label_number);
        free(essdata[i].minidisvector);
        sum_of_lab+=essdata[i].sum_of_lab;
    }



    pthread_t readthread[maxquerythread*maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter=0;
    float bsfdistance=(approximate_result.distance);

        
    readpointer.ts=ts;
    readpointer.index=index;
    readpointer.counter=&readcounter;
    readpointer.load_point=label_number;

    readpointer.lock_bsf=&lock_bsf;
    readpointer.bsf2=&bsfdistance;
    readpointer.minidisvector=minidisvector;
    readpointer.sum_of_lab=sum_of_lab;
        COUNT_CAL_TIME_END
    COUNT_INPUT_TIME_START
    for (i = 0; i < maxquerythread*maxreadthread; i++)
    {   
        pthread_create(&(readthread[i]),NULL,read_worker,(void*)&(readpointer));
    }
    for (i = 0; i < maxquerythread*maxreadthread; i++)
    {   
        pthread_join(readthread[i],NULL);
    }
    COUNT_INPUT_TIME_END
    approximate_result.distance=bsfdistance;
    free(essdata);
    free(minidisvector);
    free(label_number);
    return approximate_result;
}

pqueue_bsf exact_topk_serial_ParIS(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves, int k) 
{

    RESET_BYTES_ACCESSED

    pthread_t threadid[maxquerythread];
    pqueue_bsf *pq_bsf= pqueue_bsf_init(k);
    approximate_topk(ts, paa, index,pq_bsf);

    int sum_of_lab=0;
    // Early termination...
    if (pq_bsf->knn[k-1] == 0) {
        return *pq_bsf;
    }
    
    if(pq_bsf->knn[k-1] == FLT_MAX  || min_checked_leaves > 1) {
        refine_topk_answer(ts, paa, index, pq_bsf, minimum_distance, min_checked_leaves);
    }
    
    int i;

    SET_APPROXIMATE(pq_bsf->knn[k-1]);
    ParIS_LDCW_data *essdata=malloc(sizeof(ParIS_LDCW_data)*(maxquerythread));
    pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;
    for (i = 0; i < (maxquerythread-1); i++)
    {
        essdata[i].index=index;
        essdata[i].lock_bsf=&lock_bsf;
        essdata[i].start_number=i*(index->sax_cache_size/maxquerythread);
        essdata[i].stop_number=(i+1)*(index->sax_cache_size/maxquerythread);
        essdata[i].paa=paa;
        essdata[i].ts=ts;
        essdata[i].bsfdistance=pq_bsf->knn[k-1];
        essdata[i].sum_of_lab=0;
    }
    essdata[maxquerythread-1].index=index;
    essdata[maxquerythread-1].lock_bsf=&lock_bsf;
    essdata[maxquerythread-1].start_number=(maxquerythread-1)*(index->sax_cache_size/maxquerythread);
    essdata[maxquerythread-1].stop_number=index->sax_cache_size;
    essdata[maxquerythread-1].paa=paa;
    essdata[maxquerythread-1].ts=ts;
    essdata[maxquerythread-1].bsfdistance=pq_bsf->knn[k-1];
    essdata[maxquerythread-1].sum_of_lab=0;

    for(i=0; i<maxquerythread; i++) 
    {
        pthread_create(&(threadid[i]),NULL,mindistance_worker,(void*)&(essdata[i]));
    }
    for (i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
        sum_of_lab+=essdata[i].sum_of_lab;
    }
    unsigned long* label_number=malloc(sizeof(unsigned long)*(sum_of_lab));
    float* minidisvector=malloc(sizeof(float)*(sum_of_lab));
   
    sum_of_lab=0;

    for (i = 0; i < maxquerythread; i++)
    {
        memcpy(&(label_number[sum_of_lab]),essdata[i].label_number,sizeof(unsigned long)*essdata[i].sum_of_lab);
        memcpy(&(minidisvector[sum_of_lab]),essdata[i].minidisvector,sizeof(float)*essdata[i].sum_of_lab);
        free(essdata[i].label_number);
        free(essdata[i].minidisvector);
        sum_of_lab+=essdata[i].sum_of_lab;
    }



    pthread_t readthread[maxquerythread*maxreadthread];
    ParIS_read_worker_data readpointer;
    unsigned long readcounter=0;

        
    readpointer.ts=ts;
    readpointer.index=index;
    readpointer.counter=&readcounter;
    readpointer.load_point=label_number;

    readpointer.lock_bsf=&lock_bsf;
    readpointer.minidisvector=minidisvector;
    readpointer.sum_of_lab=sum_of_lab;
    readpointer.pq_bsf=pq_bsf;
    for (i = 0; i < maxquerythread*maxreadthread; i++)
    {   
        pthread_create(&(readthread[i]),NULL,topk_read_worker,(void*)&(readpointer));
    }
    for (i = 0; i < maxquerythread*maxreadthread; i++)
    {   
        pthread_join(readthread[i],NULL);
    }
    free(essdata);
    free(minidisvector);
    free(label_number);
    return *pq_bsf;
}

static void *mindistance_worker(void *essdata)
{    
    
    unsigned long i;
    float mindist;
    isax_index *index=((ParIS_LDCW_data*)essdata)->index;
    unsigned long start_number=((ParIS_LDCW_data*)essdata)->start_number;
    unsigned long stop_number=((ParIS_LDCW_data*)essdata)->stop_number;
    ts_type *paa=((ParIS_LDCW_data*)essdata)->paa;
    ((ParIS_LDCW_data*)essdata)->label_number=malloc(sizeof(unsigned long)*10000);
    ((ParIS_LDCW_data*)essdata)->minidisvector=malloc(sizeof(float)*10000);

    int max_number=10000;

    for(i=start_number;i<stop_number;i++)
    {

        sax_type *sax = &index->sax_cache[i * index->settings->n_segments];

        mindist = minidist_paa_to_isax_rawa_SIMD(paa, sax, index->settings->max_sax_cardinalities, index->settings);
        if(mindist <= ((ParIS_LDCW_data*)essdata)->bsfdistance) {
            if ( ((ParIS_LDCW_data*)essdata)->sum_of_lab>=max_number)
            {
                max_number += 10000;
                ((ParIS_LDCW_data*)essdata)->label_number = realloc(
                    ((ParIS_LDCW_data*)essdata)->label_number,
                    sizeof(unsigned long) * max_number);
                ((ParIS_LDCW_data*)essdata)->minidisvector = realloc(
                    ((ParIS_LDCW_data*)essdata)->minidisvector,
                    sizeof(float) * max_number);
            }
            ((ParIS_LDCW_data*)essdata)->label_number[((ParIS_LDCW_data*)essdata)->sum_of_lab]=i;
            ((ParIS_LDCW_data*)essdata)->minidisvector[((ParIS_LDCW_data*)essdata)->sum_of_lab]=mindist;
            ((ParIS_LDCW_data*)essdata)->sum_of_lab++;
        }
    }
    return NULL;
}


static void *read_worker(void *read_pointer)
{
    isax_index *index=((ParIS_read_worker_data*)read_pointer)->index;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = malloc(index->settings->ts_byte_size);
    ts_type *ts =((ParIS_read_worker_data*)read_pointer)->ts;
    unsigned long t=0,p;
    unsigned long sum_of_lab=((ParIS_read_worker_data*)read_pointer)->sum_of_lab;
    float *minidisvector=((ParIS_read_worker_data*)read_pointer)->minidisvector;
    float bsf,dist;
    while(1)
    { 
        
        pthread_rwlock_rdlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf); 
        bsf= *(((ParIS_read_worker_data*)read_pointer)->bsf2); 
        pthread_rwlock_unlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf); 
        t=__sync_fetch_and_add(((ParIS_read_worker_data*)read_pointer)->counter,1);
        if (t>=sum_of_lab) 
        {    
            break; 
        } 
        
         p=((ParIS_read_worker_data*)read_pointer)->load_point[t];
        if (minidisvector[t]<bsf) 
        {
            
            fseek(raw_file, p * index->settings->ts_byte_size, SEEK_SET); 
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
             dist = ts_euclidean_distance_SIMD(ts, ts_buffer, index->settings->timeseries_size, bsf); 
            if(dist < bsf)  
            {  
                pthread_rwlock_wrlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf); 
                if (dist<*(((ParIS_read_worker_data*)read_pointer)->bsf2)) 
                {    
                    *(((ParIS_read_worker_data*)read_pointer)->bsf2)= dist; 
                } 
                pthread_rwlock_unlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf); 
            } 
        } 
    }

    free(ts_buffer);
    fclose(raw_file);
    return NULL;
} 
static void *topk_read_worker(void *read_pointer)
{
    isax_index *index=((ParIS_read_worker_data*)read_pointer)->index;
    pqueue_bsf *pq_bsf=((ParIS_read_worker_data*)read_pointer)->pq_bsf;
    FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = malloc(index->settings->ts_byte_size);
    ts_type *ts =((ParIS_read_worker_data*)read_pointer)->ts;
    unsigned long t=0,p;
    unsigned long sum_of_lab=((ParIS_read_worker_data*)read_pointer)->sum_of_lab;
    float *minidisvector=((ParIS_read_worker_data*)read_pointer)->minidisvector;
    float bsf,dist;
    while(1)
    { 
         
        bsf= pq_bsf->knn[pq_bsf->k-1]; 
        t=__sync_fetch_and_add(((ParIS_read_worker_data*)read_pointer)->counter,1);
        if (t>=sum_of_lab) 
        {    
            break; 
        } 
        
         p=((ParIS_read_worker_data*)read_pointer)->load_point[t];
        if (minidisvector[t]<bsf) 
        {
            fseek(raw_file, p * index->settings->ts_byte_size, SEEK_SET); 
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
             dist = ts_euclidean_distance_SIMD(ts, ts_buffer, index->settings->timeseries_size, bsf); 
            if(dist <= bsf)  
            {  
                pthread_rwlock_wrlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf); 
                pqueue_bsf_insert(pq_bsf,dist,p,NULL);
                pthread_rwlock_unlock(((ParIS_read_worker_data*)read_pointer)->lock_bsf); 
            } 
        } 
    }

    free(ts_buffer);
    fclose(raw_file);
    return NULL;
}

query_result exact_search_serial_ParIS_nb(ts_type *ts, ts_type *paa, isax_index *index, float minimum_distance, int min_checked_leaves) 
{

    RESET_BYTES_ACCESSED
    pthread_t threadid[maxquerythread];
    query_result approximate_result = approximate_search(ts, paa, index);
    // Early termination...
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    if(approximate_result.distance == FLT_MAX || min_checked_leaves > 1)
    {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }
    int i;

    SET_APPROXIMATE(approximate_result.distance);
    ParIS_LDCW_data *essdata=malloc(sizeof(ParIS_LDCW_data)*(maxquerythread));
    pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;
    float bsfdistance=approximate_result.distance;
    for (i = 0; i < (maxquerythread-1); i++)
    {
        essdata[i].index=index;
        essdata[i].lock_bsf=&lock_bsf;
        essdata[i].start_number=i*(index->sax_cache_size/maxquerythread);
        essdata[i].stop_number=(i+1)*(index->sax_cache_size/maxquerythread);
        essdata[i].paa=paa;
        essdata[i].ts=ts;
        essdata[i].bsfdistance=approximate_result.distance;
        essdata[i].sum_of_lab=0;
        essdata[i].minidisvector=&bsfdistance;
    }
    essdata[maxquerythread-1].index=index;
    essdata[maxquerythread-1].lock_bsf=&lock_bsf;
    essdata[maxquerythread-1].start_number=(maxquerythread-1)*(index->sax_cache_size/maxquerythread);
    essdata[maxquerythread-1].stop_number=index->sax_cache_size;
    essdata[maxquerythread-1].paa=paa;
    essdata[maxquerythread-1].ts=ts;
    essdata[maxquerythread-1].bsfdistance=approximate_result.distance;
    essdata[maxquerythread-1].sum_of_lab=0;
    essdata[maxquerythread-1].minidisvector=&bsfdistance;

    for(i=0; i<maxquerythread; i++) 
    {
        pthread_create(&(threadid[i]),NULL,ParIS_nb_worker,(void*)&(essdata[i]));
    }
    for (i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
        if(essdata[i].bsfdistance<approximate_result.distance)
            approximate_result.distance=essdata[i].bsfdistance;
    }
    free(essdata);
    return approximate_result;
}

static void *ParIS_nb_worker(void *essdata)
{    
    unsigned long i;
    float mindist;

    isax_index *index=((ParIS_LDCW_data*)essdata)->index;
        FILE *raw_file = fopen(index->settings->raw_filename, "rb");
    fseek(raw_file, 0, SEEK_SET);
    ts_type *ts_buffer = malloc(index->settings->ts_byte_size);
    unsigned long start_number=((ParIS_LDCW_data*)essdata)->start_number;
    unsigned long stop_number=((ParIS_LDCW_data*)essdata)->stop_number;
    ts_type *paa=((ParIS_LDCW_data*)essdata)->paa;
    ts_type *ts=((ParIS_LDCW_data*)essdata)->ts;

    for(i=start_number;i<stop_number;i++)
    {

        sax_type *sax = &index->sax_cache[i * index->settings->n_segments];

        mindist = minidist_paa_to_isax_raw_SIMD(paa, sax, index->settings->max_sax_cardinalities,
                                                         index->settings);

        if(mindist <= ((ParIS_LDCW_data*)essdata)->bsfdistance) {

            fseek(raw_file, i * index->settings->ts_byte_size, SEEK_SET);
            fread(ts_buffer, index->settings->ts_byte_size, 1, raw_file);
            
            float dist = ts_euclidean_distance_SIMD(ts, ts_buffer, index->settings->timeseries_size, ((ParIS_LDCW_data*)essdata)->bsfdistance);
            if(dist<(((ParIS_LDCW_data*)essdata)->bsfdistance))
            {
                (((ParIS_LDCW_data*)essdata)->bsfdistance)=dist;
            }


        }
    }
    fclose(raw_file);
    free(ts_buffer);
    return NULL;
}


query_result exact_search_m (ts_type *ts, ts_type *paa, isax_index *index,
                           float minimum_distance, int min_checked_leaves) 
{
    query_result approximate_result = approximate_search_SIMD(ts, paa, index);
    query_result bsf_result = approximate_result;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    int i;
    if (approximate_result.distance == 0) {
        return approximate_result;
    }
    if(approximate_result.distance == FLT_MAX || min_checked_leaves > 1) {
        approximate_result = refine_answer(ts, paa, index, approximate_result, minimum_distance, min_checked_leaves);
    }

    pqueue_t *pq = pqueue_init(index->settings->root_nodes_size,
                               cmp_pri, get_pri, set_pri, get_pos, set_pos);



    query_result *do_not_remove = &approximate_result;

    SET_APPROXIMATE(approximate_result.distance);

    RESET_BYTES_ACCESSED

    if(approximate_result.node != NULL) {
        pqueue_insert(pq, &approximate_result);
    }

    pthread_t threadid[maxquerythread];

    refind_answer_fonction_data rfdata;
    pthread_mutex_t lock_queue=PTHREAD_MUTEX_INITIALIZER,lock_current_root_node=PTHREAD_MUTEX_INITIALIZER;
    pthread_rwlock_t lock_bsf=PTHREAD_RWLOCK_INITIALIZER;
    pthread_barrier_t lock_barrier;
    pthread_barrier_init(&lock_barrier, NULL, maxquerythread);
    rfdata.paa=paa;
    rfdata.ts=ts;
    
    rfdata.lock_queue=&lock_queue;
    rfdata.lock_current_root_node=&lock_current_root_node;
    rfdata.lock_bsf=&lock_bsf;
    
    rfdata.index=index;
    rfdata.minimum_distance=minimum_distance;
    rfdata.pq=pq;
    rfdata.current_root_node = index->first_node;
    rfdata.bsf_result = &bsf_result;
   
    query_result * n;
    
    rfdata.lock_barrier=&lock_barrier;
    for (i = 0; i < maxquerythread; i++)
    {
        pthread_create(&(threadid[i]),NULL,exact_search_fonction,(void*)&(rfdata));
    }
    for (i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
    }


    while ((n = pqueue_pop(pq)))
    {
        if(n != do_not_remove)
            free(n);
    }
    pthread_barrier_destroy(&lock_barrier);
    pqueue_free(pq);
    
    return *(rfdata.bsf_result);
}
static void *exact_search_fonction(void *rfdata)
{
    isax_node *current_root_node;
    query_result *n;
    isax_index *index=((refind_answer_fonction_data*)rfdata)->index;
    ts_type *paa=((refind_answer_fonction_data*)rfdata)->paa;
    ts_type *ts=((refind_answer_fonction_data*)rfdata)->ts;
    pqueue_t *pq=((refind_answer_fonction_data*)rfdata)->pq;
    query_result *do_not_remove = ((refind_answer_fonction_data*)rfdata)->bsf_result;
    float minimum_distance=((refind_answer_fonction_data*)rfdata)->minimum_distance;
    float bsfdisntance;
    int tight_bound = index->settings->tight_bound;
    int aggressive_check = index->settings->aggressive_check;
    query_result *bsf_result=(((refind_answer_fonction_data*)rfdata)->bsf_result);
    while (1) 
    {
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
        current_root_node= ((refind_answer_fonction_data*)rfdata)->current_root_node;
        if (current_root_node != NULL)
        {
            ((refind_answer_fonction_data*)rfdata)->current_root_node=((refind_answer_fonction_data*)rfdata)->current_root_node->next;
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
        }
        else    
        {
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_current_root_node);
            break;
        }
        
        

        query_result * mindist_result = malloc(sizeof(query_result));
        mindist_result->distance =  minidist_paa_to_isax(paa, current_root_node->isax_values, current_root_node->isax_cardinalities, index->settings, 0);

        mindist_result->node = current_root_node;
        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);

        pqueue_insert(pq, mindist_result);

        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
    }
    pthread_barrier_wait(((refind_answer_fonction_data*)rfdata)->lock_barrier);

 while (1)
    {


        pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
        n = pqueue_pop(pq);
        pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
        if(n==NULL)
            break;
        pthread_rwlock_rdlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
        bsfdisntance=bsf_result->distance;
        pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
        // The best node has a worse mindist, so search is finished!
        if (n->distance >= bsfdisntance || n->distance > minimum_distance) {
            pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
            pqueue_insert(pq, n);
            pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
            break;
        }
        else 
        {
            // If it is a leaf, check its real distance.
            if (n->node->is_leaf) {
                if (!n->node->has_full_data_file &&
                    (n->node->leaf_size > index->settings->min_leaf_size))
                {
                    split_node(index, n->node, 0, 1);
                    pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    pqueue_insert(pq, n);
                    pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    continue;
                }
                if(tight_bound) {
                    float mindistance = calculate_minimum_distance_SIMD(index, n->node, ts, paa);
                    if(mindistance >= bsfdisntance)
                    {
                        if(n != do_not_remove)
                            free(n);
                        continue;
                    }
                }
                float distance = calculate_node_distance_SIMD(index, n->node, ts, bsfdisntance);
                if (distance < bsfdisntance)
                {
                    pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                    if(distance < bsf_result->distance)
                    {
                        bsf_result->distance = distance;
                        bsf_result->node = n->node;
                    }

                    pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                }
            }
            else {
                // If it is an intermediate node calculate mindist for children
                // and push them in the queue
                if (n->node->left_child->isax_cardinalities != NULL) {
                    if(n->node->left_child->is_leaf && !n->node->left_child->has_partial_data_file && aggressive_check){
                        float distance = calculate_node_distance_SIMD(index, n->node->left_child, ts, bsfdisntance);
                        if (distance < bsfdisntance)
                        {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                            if(distance <bsf_result->distance)
                            {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->left_child;
                            }
                            pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                        }
                    }
                    else {
                    query_result * mindist_result = malloc(sizeof(query_result));
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->left_child->isax_values, n->node->left_child->isax_cardinalities, index->settings, 0);
                    mindist_result->node = n->node->left_child;
                    pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    pqueue_insert(pq, mindist_result);
                    pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    }
                }
                if (n->node->right_child->isax_cardinalities != NULL) {
                    if(n->node->right_child->is_leaf && !n->node->right_child->has_partial_data_file && aggressive_check){
                        float distance = calculate_node_distance_SIMD(index, n->node->right_child, ts, bsfdisntance);
                        if (distance < bsfdisntance)
                        {
                            pthread_rwlock_wrlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                            if(distance <bsf_result->distance)
                            {
                                bsf_result->distance = distance;
                                bsf_result->node = n->node->right_child;
                            }
                            pthread_rwlock_unlock(((refind_answer_fonction_data*)rfdata)->lock_bsf);
                        }
                    }
                    else {
                    query_result * mindist_result = malloc(sizeof(query_result));
                    mindist_result->distance =  minidist_paa_to_isax(paa, n->node->right_child->isax_values, n->node->right_child->isax_cardinalities, index->settings, 0);
                    mindist_result->node = n->node->right_child;
                    pthread_mutex_lock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    pqueue_insert(pq, mindist_result);
                    pthread_mutex_unlock(((refind_answer_fonction_data*)rfdata)->lock_queue);
                    }
                }
            }

            // Free the node currently popped.
            if(n != do_not_remove)
                free(n);
        }
    }
    return NULL;
}
