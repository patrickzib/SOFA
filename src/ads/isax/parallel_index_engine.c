//
//  parallel_index_engine.c
//
//
//  Created by Botao PENG on 29/1/18.
//


#include "config.h"
#include "../../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <pthread.h>
#include <unistd.h>
#include <stdbool.h>

#include "ads/isax_node.h"
#include "ads/isax_index.h"
#include "ads/isax_query_engine.h"
#include "ads/isax_node_record.h"
#include "ads/isax_file_loaders.h"
#include "ads/isax_first_buffer_layer.h"
#include "ads/parallel_index_engine.h"
#include "ads/parallel_query_engine.h"

int read_block_length;

typedef struct {
    file_position_type pos;
    isax_index *index;
    ts_type *ts;
    sax_type *saxv;
    int fin_number;
    int blocid;
    pthread_mutex_t *lock_cbl;
    pthread_mutex_t *lock_firstnode;
    int *bufferpresize;
    int *nodecounter;
    pthread_barrier_t *lock_barrier1;
    pthread_barrier_t *lock_barrier2;
    pthread_barrier_t *lock_barrier3;
    bool finished;
} index_buffer_data;

static void *indexbulkloadingworker(void *transferdata);
static void *indexbulkloadingworker_new(void *transferdata);
static isax_node *insert_to_fbl_m(first_buffer_layer *fbl, sax_type *sax,
                                  file_position_type *pos, root_mask_type mask,
                                  isax_index *index, pthread_mutex_t *lock_firstnode);
static enum response indexconstruction(first_buffer_layer *fbl, isax_index *index,
                                       pthread_mutex_t *lock_index, int calculate_thread);
static enum response indexflush(first_buffer_layer *fbl, isax_index *index,
                                pthread_mutex_t *lock_index, int calculate_thread);
static void *indexconstructionworker(void *input);
static void *indexflushworker(void *input);
static root_mask_type isax_fbl_index_insert_m(isax_index *index, sax_type *sax,
                                              file_position_type *pos,
                                              pthread_mutex_t *lock_cbl,
                                              pthread_mutex_t *lock_firstnode);


void isax_index_binary_file_m(const char *ifilename, int ts_num, isax_index *index,int calculate_thread)
{
    fprintf(stderr, ">>> Indexing: %s\n", ifilename);
    FILE * ifile;

    COUNT_INPUT_TIME_START
    ifile = fopen (ifilename,"rb");
    COUNT_INPUT_TIME_END
    index_buffer_data *input_data=malloc(sizeof(index_buffer_data)*(calculate_thread-1));
    file_position_type *pos = malloc(sizeof(file_position_type));
    pthread_t threadid[calculate_thread-1];
    bool sax_fist_time_check = false;
    int j,conter=0;
    long long int i;
    int prev_flush_time=0,now_flush_time=0;
    int sax_save_number;
    //initial the locks
    pthread_mutex_t lock_index=PTHREAD_MUTEX_INITIALIZER,
                    lock_firstnode=PTHREAD_MUTEX_INITIALIZER,lock_disk=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t *lockcbl;
    pthread_barrier_t lock_barrier1, lock_barrier2;
    pthread_barrier_init(&lock_barrier1, NULL, calculate_thread);
    pthread_barrier_init(&lock_barrier2, NULL, calculate_thread);

    // set the thread on decided cpu

    lockcbl= malloc(sizeof(pthread_mutex_t)*LOCK_SIZE);

    for(i=0; i<LOCK_SIZE ;i++)
    {
        pthread_mutex_init(&lockcbl[i], NULL);
    }

    for (i = 0; i < (calculate_thread-1); i++)
    {
        input_data[i].index=index;
        input_data[i].lock_cbl=lockcbl;
        input_data[i].lock_firstnode =&lock_firstnode;
        input_data[i].blocid=i*read_block_length;

        input_data[i].lock_barrier1=&lock_barrier1;
        input_data[i].lock_barrier2=&lock_barrier2;
        input_data[i].finished=false;
    }
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }

    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz/index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);

    if (ts_num < 0 || total_records < (file_position_type)ts_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }


    ts_type * ts     = malloc(sizeof(ts_type) * index->settings->timeseries_size*read_block_length*(calculate_thread-1));
    ts_type * ts1    = malloc(sizeof(ts_type) * index->settings->timeseries_size*read_block_length*(calculate_thread-1));
    ts_type * ts2;
    sax_type * saxv  = malloc(sizeof(sax_type) * index->settings->n_segments*read_block_length*(calculate_thread-1));
    sax_type * saxv1 = malloc(sizeof(sax_type) * index->settings->n_segments*read_block_length*(calculate_thread-1));
    sax_type * saxv2;

    index->settings->raw_filename = malloc(256);
    strcpy(index->settings->raw_filename, ifilename);

    *pos = ftell(ifile);


    if(ts_num>read_block_length*(calculate_thread-1))
    {   COUNT_INPUT_TIME_START
        fread(ts1, sizeof(ts_type), index->settings->timeseries_size*read_block_length*(calculate_thread-1), ifile);
        COUNT_INPUT_TIME_END
        ts2   =ts;
        ts    =ts1;
        ts1   =ts2;

        for ( j = 0; j < (calculate_thread-1); j++)
        {
            input_data[j].pos        =*pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
            input_data[j].ts         =&(ts[index->settings->timeseries_size*j*read_block_length]);
            input_data[j].saxv       =saxv;
            input_data[j].fin_number =read_block_length;
        }
        for (j = 0; j < (calculate_thread-1); j++)
        {
            pthread_create(&(threadid[j]),NULL,indexbulkloadingworker,(void*)&(input_data[j]));
        }

        for (i = read_block_length*(calculate_thread-1)*2; i <= ts_num; i+=read_block_length*(calculate_thread-1))
        {

            *pos = ftell(ifile);
            //read the data of next round
            COUNT_INPUT_TIME_START
            fread(ts1, sizeof(ts_type), index->settings->timeseries_size*read_block_length*(calculate_thread-1), ifile);
            COUNT_INPUT_TIME_END
            //write the sax in disk (last round)
            if (sax_fist_time_check)
            {
                COUNT_OUTPUT_TIME_START
                fwrite(saxv1, index->settings->sax_byte_size, read_block_length*(calculate_thread-1), index->sax_file);
                COUNT_OUTPUT_TIME_END
            }
            else
            {
                sax_fist_time_check=true;
            }

            pthread_barrier_wait(&lock_barrier1);
            //wait for the finish of other threads
            ts2   =ts;
            ts    =ts1;
            ts1   =ts2;
            __sync_fetch_and_add(&(index->fbl->current_record_index),read_block_length*(calculate_thread-1));
            now_flush_time=i/(index->settings->max_total_buffer_size);
            if(now_flush_time!=prev_flush_time)
            {
                COUNT_QUEUE_TIME_START
                indexconstruction(index->fbl, index, &lock_index, calculate_thread);
                COUNT_QUEUE_TIME_END
            }
            saxv2 =saxv;
            saxv  =saxv1;
            saxv1 =saxv2;

            prev_flush_time=now_flush_time;
            for ( j = 0; j < (calculate_thread-1); j++)
            {
                input_data[j].pos        =*pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
                input_data[j].ts         =&(ts[index->settings->timeseries_size*j*read_block_length]);
                input_data[j].saxv       =saxv;
                input_data[j].fin_number =read_block_length;
            }
            pthread_barrier_wait(&lock_barrier2);
        }

        pthread_barrier_wait(&lock_barrier1);
        for (j = 0; j < (calculate_thread-1); j++)
        {
            input_data[j].finished=true;
        }
            //wait for the finish of other threads
        __sync_fetch_and_add(&(index->fbl->current_record_index),read_block_length*(calculate_thread-1));
        now_flush_time=i/(index->settings->max_total_buffer_size);
        if(now_flush_time!=prev_flush_time)
        {
            COUNT_QUEUE_TIME_START
            indexconstruction(index->fbl, index, &lock_index, calculate_thread);
            COUNT_QUEUE_TIME_END
        }

        prev_flush_time=now_flush_time;
        for ( j = 0; j < (calculate_thread-1); j++)
        {
            input_data[j].pos        =*pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
            input_data[j].ts         =&(ts[index->settings->timeseries_size*j*read_block_length]);
            input_data[j].saxv       =saxv;
            input_data[j].fin_number =read_block_length;
        }
        pthread_barrier_wait(&lock_barrier2);
        for (j = 0; j < (calculate_thread-1); j++)
        {
            pthread_join(threadid[j],NULL);
        }
    }
    COUNT_INPUT_TIME_START
    *pos = ftell(ifile);
    fread(ts1, sizeof(ts_type), index->settings->timeseries_size*(ts_num%(read_block_length*(calculate_thread-1))), ifile);
    COUNT_INPUT_TIME_END
    ts2   =ts;
    ts    =ts1;
    ts1   =ts2;



    //handle the rest data
    int conter_ts_number=ts_num%(read_block_length*(calculate_thread-1));
    sax_save_number=conter_ts_number;

    for ( j = 0; j < (calculate_thread-1); j++)
    {
        input_data[j].pos = *pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
        conter++;
        input_data[j].ts=&(ts[index->settings->timeseries_size*j*read_block_length]);
        input_data[j].fin_number=min(conter_ts_number,read_block_length);
        input_data[j].saxv=saxv;
        conter_ts_number=conter_ts_number-read_block_length;

        if(conter_ts_number<0)
            break;
    }

    pthread_barrier_init(&lock_barrier1, NULL, conter+1);
    pthread_barrier_init(&lock_barrier2, NULL, conter+1);
    for ( j = 0; j < conter; j++)
    {
        input_data[j].finished=false;
        pthread_create(&(threadid[j]),NULL,indexbulkloadingworker,(void*)&(input_data[j]));
    }

    if (sax_fist_time_check)
    {
        pthread_mutex_lock(&lock_disk);
        COUNT_OUTPUT_TIME_START
        fwrite(saxv1, index->settings->sax_byte_size, read_block_length*(calculate_thread-1), index->sax_file);
        COUNT_OUTPUT_TIME_END
        pthread_mutex_unlock(&lock_disk);
    }
    else
    {
        sax_fist_time_check=true;
    }
    pthread_barrier_wait(&lock_barrier1);
    saxv2 =saxv;
    saxv  =saxv1;
    saxv1 =saxv2;
    for ( j = 0; j < conter; j++)
    {
        input_data[j].finished=true;
    }
    pthread_barrier_wait(&lock_barrier2);

    for (j = 0; j < conter; j++)
    {
        pthread_join(threadid[j],NULL);
    }
    pthread_mutex_lock(&lock_disk);
    COUNT_OUTPUT_TIME_START
    fwrite(saxv1, index->settings->sax_byte_size, sax_save_number, index->sax_file);
    COUNT_OUTPUT_TIME_END
    pthread_mutex_unlock(&lock_disk);
    __sync_fetch_and_add(&(index->fbl->current_record_index),sax_save_number);
    COUNT_QUEUE_TIME_START
    indexconstruction(index->fbl, index, &lock_index, calculate_thread);
    COUNT_QUEUE_TIME_END
    free(ts);
    free(ts1);
    free(input_data);
    free(lockcbl);
    free(saxv);
    free(saxv1);
    free(pos);
    pthread_barrier_destroy(&lock_barrier1);
    pthread_barrier_destroy(&lock_barrier2);
    COUNT_INPUT_TIME_START
    fclose(ifile);
    COUNT_INPUT_TIME_END
}

void isax_index_binary_file_m_new(const char *ifilename, int ts_num, isax_index *index,int calculate_thread)
{
    fprintf(stderr, ">>> Indexing: %s\n", ifilename);
    FILE * ifile;
    COUNT_INPUT_TIME_START
    ifile = fopen (ifilename,"rb");
    COUNT_INPUT_TIME_END
    index_buffer_data *input_data=malloc(sizeof(index_buffer_data)*(calculate_thread-1));
    file_position_type *pos = malloc(sizeof(file_position_type));
    pthread_t threadid[calculate_thread-1];
    bool sax_fist_time_check = false;
    int j,conter=0;
    long long int i;
    int prev_flush_time=0,now_flush_time=0;
    int sax_save_number;
    int nodecounter=0;
    int bufferpresize[index->fbl->number_of_buffers];

    for (int i = 0; i < index->fbl->number_of_buffers; i++)
    {
        bufferpresize[i]=0;
    }
    //initial the locks
    pthread_mutex_t lock_index=PTHREAD_MUTEX_INITIALIZER,
                    lock_firstnode=PTHREAD_MUTEX_INITIALIZER,lock_disk=PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t *lockcbl;
    pthread_barrier_t lock_barrier1, lock_barrier2, lock_barrier3;
    pthread_barrier_init(&lock_barrier1, NULL, calculate_thread);
    pthread_barrier_init(&lock_barrier2, NULL, calculate_thread);
    pthread_barrier_init(&lock_barrier3, NULL, calculate_thread-1);
    // set the thread on decided cpu

    lockcbl= malloc(sizeof(pthread_mutex_t)*LOCK_SIZE);

    for(i=0; i<LOCK_SIZE ;i++)
    {
        pthread_mutex_init(&lockcbl[i], NULL);
    }

    for (i = 0; i < (calculate_thread-1); i++)
    {
        input_data[i].index=index;
        input_data[i].lock_cbl=lockcbl;
        input_data[i].lock_firstnode =&lock_firstnode;
        input_data[i].blocid=i*read_block_length;
        input_data[i].bufferpresize=bufferpresize;
        input_data[i].lock_barrier1=&lock_barrier1;
        input_data[i].lock_barrier2=&lock_barrier2;
        input_data[i].lock_barrier3=&lock_barrier3;
        input_data[i].finished=false;
        input_data[i].nodecounter=&nodecounter;
    }
    if (ifile == NULL) {
        fprintf(stderr, "File %s not found!\n", ifilename);
        exit(-1);
    }

    fseek(ifile, 0L, SEEK_END);
    file_position_type sz = (file_position_type) ftell(ifile);
    file_position_type total_records = sz/index->settings->ts_byte_size;
    fseek(ifile, 0L, SEEK_SET);

    if (ts_num < 0 || total_records < (file_position_type)ts_num) {
        fprintf(stderr, "File %s has only %llu records!\n", ifilename, total_records);
        exit(-1);
    }


    ts_type * ts     = malloc(sizeof(ts_type) * index->settings->timeseries_size*read_block_length*(calculate_thread-1));
    ts_type * ts1    = malloc(sizeof(ts_type) * index->settings->timeseries_size*read_block_length*(calculate_thread-1));
    ts_type * ts2;
    sax_type * saxv  = malloc(sizeof(sax_type) * index->settings->n_segments*read_block_length*(calculate_thread-1));
    sax_type * saxv1 = malloc(sizeof(sax_type) * index->settings->n_segments*read_block_length*(calculate_thread-1));
    sax_type * saxv2;

    index->settings->raw_filename = malloc(256);
    strcpy(index->settings->raw_filename, ifilename);


    *pos = ftell(ifile);


    if(ts_num>read_block_length*(calculate_thread-1))
    {
        COUNT_INPUT_TIME_START
        fread(ts1, sizeof(ts_type), index->settings->timeseries_size*read_block_length*(calculate_thread-1), ifile);
        COUNT_INPUT_TIME_END
        ts2   =ts;
        ts    =ts1;
        ts1   =ts2;

        for ( j = 0; j < (calculate_thread-1); j++)
        {
            input_data[j].pos        =*pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
            input_data[j].ts         =&(ts[index->settings->timeseries_size*j*read_block_length]);
            input_data[j].saxv       =saxv;
            input_data[j].fin_number =read_block_length;
        }
        for (j = 0; j < (calculate_thread-1); j++)
        {
            pthread_create(&(threadid[j]),NULL,indexbulkloadingworker_new,(void*)&(input_data[j]));
        }

        for (i = read_block_length*(calculate_thread-1)*2; i <= ts_num; i+=read_block_length*(calculate_thread-1))
        {

            *pos = ftell(ifile);
            //read the data of next round
            COUNT_INPUT_TIME_START
            fread(ts1, sizeof(ts_type), index->settings->timeseries_size*read_block_length*(calculate_thread-1), ifile);
            COUNT_INPUT_TIME_END

            //write the sax in disk (last round)
            if (sax_fist_time_check)
            {
                COUNT_OUTPUT_TIME_START
                fwrite(saxv1, index->settings->sax_byte_size, read_block_length*(calculate_thread-1), index->sax_file);
                COUNT_OUTPUT_TIME_END
            }
            else
            {
                sax_fist_time_check=true;
            }

            pthread_barrier_wait(&lock_barrier1);
            //wait for the finish of other threads
            ts2   =ts;
            ts    =ts1;
            ts1   =ts2;
            __sync_fetch_and_add(&(index->fbl->current_record_index),read_block_length*(calculate_thread-1));
            now_flush_time=i/(index->settings->max_total_buffer_size);
            if(now_flush_time!=prev_flush_time)
            {
                COUNT_OUTPUT_TIME_START
                indexflush(index->fbl, index, &lock_index, calculate_thread);
                COUNT_OUTPUT_TIME_END
                    for (int i = 0; i < index->fbl->number_of_buffers; i++)
            {
                    bufferpresize[i]=0;
            }
            }
            saxv2 =saxv;
            saxv  =saxv1;
            saxv1 =saxv2;

            prev_flush_time=now_flush_time;
            for ( j = 0; j < (calculate_thread-1); j++)
            {
                input_data[j].pos        =*pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
                input_data[j].ts         =&(ts[index->settings->timeseries_size*j*read_block_length]);
                input_data[j].saxv       =saxv;
                input_data[j].fin_number =read_block_length;
            }
            nodecounter=0;
            pthread_barrier_wait(&lock_barrier2);

        }

        pthread_barrier_wait(&lock_barrier1);
        for (j = 0; j < (calculate_thread-1); j++)
        {
            input_data[j].finished=true;
        }
            //wait for the finish of other threads
        __sync_fetch_and_add(&(index->fbl->current_record_index),read_block_length*(calculate_thread-1));
        now_flush_time=i/(index->settings->max_total_buffer_size);
        if(now_flush_time!=prev_flush_time)
        {
            COUNT_OUTPUT_TIME_START
            indexflush(index->fbl, index, &lock_index, calculate_thread);
            COUNT_OUTPUT_TIME_END
        for (int i = 0; i < index->fbl->number_of_buffers; i++)
            {
                bufferpresize[i]=0;
            }
        }

        prev_flush_time=now_flush_time;
        for ( j = 0; j < (calculate_thread-1); j++)
        {
            input_data[j].pos        =*pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
            input_data[j].ts         =&(ts[index->settings->timeseries_size*j*read_block_length]);
            input_data[j].saxv       =saxv;
            input_data[j].fin_number =read_block_length;
        }
        nodecounter=0;
        pthread_barrier_wait(&lock_barrier2);
        for (j = 0; j < (calculate_thread-1); j++)
        {
            pthread_join(threadid[j],NULL);
        }
    }
    *pos = ftell(ifile);
    COUNT_INPUT_TIME_START
    fread(ts1, sizeof(ts_type), index->settings->timeseries_size*(ts_num%(read_block_length*(calculate_thread-1))), ifile);
    COUNT_INPUT_TIME_END
    ts2   =ts;
    ts    =ts1;
    ts1   =ts2;



    //handle the rest data
    int conter_ts_number=ts_num%(read_block_length*(calculate_thread-1));
    sax_save_number=conter_ts_number;

    for ( j = 0; j < (calculate_thread-1); j++)
    {
        input_data[j].pos = *pos+index->settings->timeseries_size*sizeof(ts_type)*j*read_block_length;
        conter++;
        input_data[j].ts=&(ts[index->settings->timeseries_size*j*read_block_length]);
        input_data[j].fin_number=min(conter_ts_number,read_block_length);
        input_data[j].saxv=saxv;
        conter_ts_number=conter_ts_number-read_block_length;

        if(conter_ts_number<0)
            break;
    }

    pthread_barrier_init(&lock_barrier1, NULL, conter+1);
    pthread_barrier_init(&lock_barrier2, NULL, conter+1);
    pthread_barrier_init(&lock_barrier3, NULL, conter);
    for ( j = 0; j < conter; j++)
    {
        input_data[j].finished=false;
        pthread_create(&(threadid[j]),NULL,indexbulkloadingworker_new,(void*)&(input_data[j]));
    }

    if (sax_fist_time_check)
    {
        pthread_mutex_lock(&lock_disk);
        COUNT_OUTPUT_TIME_START
        fwrite(saxv1, index->settings->sax_byte_size, read_block_length*(calculate_thread-1), index->sax_file);
        COUNT_OUTPUT_TIME_END
        pthread_mutex_unlock(&lock_disk);
    }
    else
    {
        sax_fist_time_check=true;
    }

    pthread_barrier_wait(&lock_barrier1);
    saxv2 =saxv;
    saxv  =saxv1;
    saxv1 =saxv2;
    for ( j = 0; j < conter; j++)
    {
        input_data[j].finished=true;
    }
    nodecounter=0;
    pthread_barrier_wait(&lock_barrier2);

    for (j = 0; j < conter; j++)
    {
        pthread_join(threadid[j],NULL);
    }
    pthread_mutex_lock(&lock_disk);
    COUNT_OUTPUT_TIME_START
    fwrite(saxv1, index->settings->sax_byte_size, sax_save_number, index->sax_file);
    COUNT_OUTPUT_TIME_END
    pthread_mutex_unlock(&lock_disk);
    __sync_fetch_and_add(&(index->fbl->current_record_index),sax_save_number);
    COUNT_OUTPUT_TIME_START
    indexflush(index->fbl, index, &lock_index, calculate_thread);
    COUNT_OUTPUT_TIME_END
        for (int i = 0; i < index->fbl->number_of_buffers; i++)
    {
        bufferpresize[i]=0;
    }
    free(ts);
    free(ts1);
    free(input_data);
    free(lockcbl);
    free(saxv);
    free(saxv1);
    free(pos);
    pthread_barrier_destroy(&lock_barrier1);
    pthread_barrier_destroy(&lock_barrier2);
    pthread_barrier_destroy(&lock_barrier3);
    COUNT_INPUT_TIME_START
    fclose(ifile);
    COUNT_INPUT_TIME_END
}

static void *indexbulkloadingworker(void *transferdata)
{
    sax_type * sax = malloc(sizeof(sax_type) * ((index_buffer_data*)transferdata)->index->settings->n_segments);

    int fin_number=((index_buffer_data*)transferdata)->fin_number;
    file_position_type *pos = malloc(sizeof(file_position_type));
    sax_type * saxv;
    int offset_saxv=((index_buffer_data*)transferdata)->blocid*((index_buffer_data*)transferdata)->index->settings->n_segments;
    int n_segments=((index_buffer_data*)transferdata)->index->settings->n_segments;
    isax_index *index= ((index_buffer_data*)transferdata)->index;
    int i=0;
    pthread_barrier_t *lock_barrier1, *lock_barrier2;

    while(!((index_buffer_data*)transferdata)->finished)
    {
        saxv=(((index_buffer_data*)transferdata)->saxv);
        lock_barrier1=((index_buffer_data*)transferdata)->lock_barrier1;
        lock_barrier2=((index_buffer_data*)transferdata)->lock_barrier2;
        for (i=0;i<fin_number;i++)
        {
            if(sax_from_ts(((ts_type*)(((index_buffer_data*)transferdata)->ts)+i*index->settings->timeseries_size), sax, index->settings) == SUCCESS)
            {
                *pos= ((index_buffer_data*)transferdata)->pos+index->settings->timeseries_size*sizeof(ts_type)*i;
                memcpy(&(saxv[offset_saxv+i*n_segments]), sax, sizeof(sax_type) * n_segments);
                isax_fbl_index_insert_m(index, sax, pos,
                                        ((index_buffer_data*)transferdata)->lock_cbl,
                                        ((index_buffer_data*)transferdata)->lock_firstnode);
            }
            else
            {
                fprintf(stderr, "error: cannot insert record in index, since sax representation\
                    failed to be created");
            }

        }
        pthread_barrier_wait(lock_barrier1);
        pthread_barrier_wait(lock_barrier2);
    }
    free(pos);
    free(sax);
    return NULL;
}
static void *indexbulkloadingworker_new(void *transferdata)
{
    sax_type * sax = malloc(sizeof(sax_type) * ((index_buffer_data*)transferdata)->index->settings->n_segments);

    int fin_number=((index_buffer_data*)transferdata)->fin_number;
    file_position_type *pos = malloc(sizeof(file_position_type));
    sax_type * saxv;
    int offset_saxv=((index_buffer_data*)transferdata)->blocid*((index_buffer_data*)transferdata)->index->settings->n_segments;
    int n_segments=((index_buffer_data*)transferdata)->index->settings->n_segments;
    isax_index *index= ((index_buffer_data*)transferdata)->index;
    int i=0,j;
    pthread_barrier_t *lock_barrier1=((index_buffer_data*)transferdata)->lock_barrier1;
    pthread_barrier_t *lock_barrier2=((index_buffer_data*)transferdata)->lock_barrier2;
    pthread_barrier_t *lock_barrier3=((index_buffer_data*)transferdata)->lock_barrier3;
    isax_node_record *r = calloc(1, sizeof(isax_node_record));
    while(!((index_buffer_data*)transferdata)->finished)
    {
        saxv=(((index_buffer_data*)transferdata)->saxv);

        for (i=0;i<fin_number;i++)
        {
            if(sax_from_ts(((ts_type*)(((index_buffer_data*)transferdata)->ts)+i*index->settings->timeseries_size), sax, index->settings) == SUCCESS)
            {
                *pos= ((index_buffer_data*)transferdata)->pos+index->settings->timeseries_size*sizeof(ts_type)*i;
                memcpy(&(saxv[offset_saxv+i*n_segments]), sax, sizeof(sax_type) * n_segments);
                isax_fbl_index_insert_m(index, sax, pos,
                                        ((index_buffer_data*)transferdata)->lock_cbl,
                                        ((index_buffer_data*)transferdata)->lock_firstnode);
            }
            else
            {
                fprintf(stderr, "error: cannot insert record in index, since sax representation\
                    failed to be created");
            }

        }
        pthread_barrier_wait(lock_barrier3);
        while(1)
        {
            j=__sync_fetch_and_add(((index_buffer_data*)transferdata)->nodecounter,1);
            if(j>=index->fbl->number_of_buffers)
            {
                break;
            }
            fbl_soft_buffer *current_fbl_node = &index->fbl->soft_buffers[j];
            if (!current_fbl_node->initialized) {
                continue;
            }
            for (i=((index_buffer_data*)transferdata)->bufferpresize[j]; i<current_fbl_node->buffer_size; i++)
            {
                r->sax = (sax_type *) current_fbl_node->sax_records[i];
                r->position = (file_position_type *) current_fbl_node->pos_records[i];
                r->insertion_mode = NO_TMP | PARTIAL;

                // Add record to index
                add_record_to_node(index, current_fbl_node->node, r, 1, 1);

            }
            ((index_buffer_data*)transferdata)->bufferpresize[j]=current_fbl_node->buffer_size;
        }

        pthread_barrier_wait(lock_barrier1);
        pthread_barrier_wait(lock_barrier2);
    }
    free(pos);
    free(sax);
    free(r);
    return NULL;
}
static isax_node *insert_to_fbl_m(first_buffer_layer *fbl, sax_type *sax,
                          file_position_type *pos,root_mask_type mask,
                          isax_index *index, pthread_mutex_t *lock_firstnode)
{
    fbl_soft_buffer *current_buffer = &fbl->soft_buffers[(int) mask];
    char *cd_s, *cd_p;
    // Check if this buffer is initialized
    if (!current_buffer->initialized) {
    #ifdef DEBUG
        printf("*** Creating new FBL node. ***\n\n");
    #endif
        current_buffer->initialized = 1;
        current_buffer->max_buffer_size = 0;
        current_buffer->buffer_size = 0;

        current_buffer->node = isax_root_node_init(mask,
                                                   index->settings->initial_leaf_buffer_size);

        pthread_mutex_lock(lock_firstnode);

        __sync_fetch_and_add(&(index->root_nodes), 1);


        if(index->first_node == NULL)
        {
            index->first_node = current_buffer->node;
            pthread_mutex_unlock(lock_firstnode);
            current_buffer->node->next = NULL;
            current_buffer->node->previous = NULL;

        }
        else
        {
            isax_node * prev_first = index->first_node;
            index->first_node = current_buffer->node;
            index->first_node->next = prev_first;

            prev_first->previous = current_buffer->node;
            pthread_mutex_unlock(lock_firstnode);
        }
    }

    // Check if this buffer is not full!
    if (current_buffer->buffer_size >= current_buffer->max_buffer_size) {

        if(current_buffer->max_buffer_size == 0) {
            current_buffer->max_buffer_size = fbl->initial_buffer_size;

            current_buffer->sax_records = malloc(sizeof(sax_type *) *
                                                 current_buffer->max_buffer_size);
            current_buffer->pos_records = malloc(sizeof(file_position_type *) *
                                                 current_buffer->max_buffer_size);
        }
        else {
            current_buffer->max_buffer_size *= BUFFER_REALLOCATION_RATE;

            current_buffer->sax_records = realloc(current_buffer->sax_records,
                                           sizeof(sax_type *) *
                                           current_buffer->max_buffer_size);
            current_buffer->pos_records = realloc(current_buffer->pos_records,
                                           sizeof(file_position_type *) *
                                           current_buffer->max_buffer_size);

        }
    }
    if (current_buffer->sax_records == NULL || current_buffer->pos_records == NULL) {
        fprintf(stderr, "error: Could not allocate memory in FBL.");
        return OUT_OF_MEMORY_FAILURE;
    }

    // Copy data to hard buffer and make current buffer point to the hard one
    cd_s= __sync_fetch_and_add(&(fbl->current_record),index->settings->sax_byte_size);

    cd_p= __sync_fetch_and_add(&(fbl->current_record),index->settings->position_byte_size);

    current_buffer->sax_records[current_buffer->buffer_size] = (sax_type*) cd_s;
    memcpy((void *) cd_s, (void *) sax, index->settings->sax_byte_size);


    current_buffer->pos_records[current_buffer->buffer_size] = (file_position_type*) cd_p;
    memcpy((void *) cd_p, (void *) pos, index->settings->position_byte_size);





    #ifdef DEBUG
    printf("*** Added to node ***\n\n");
    #ifdef TOY
    sax_print(sax, index->settings->n_segments,
              index->settings->sax_bit_cardinality);
    #endif
    #endif
    current_buffer->buffer_size++;
    return current_buffer->node;
}
isax_node * insert_to_pRecBuf(parallel_first_buffer_layer *fbl, sax_type *sax,
                          file_position_type *pos,root_mask_type mask,
                          isax_index *index, pthread_mutex_t *lock_firstnode, int workernumber,int total_workernumber)
{
    parallel_fbl_soft_buffer *current_buffer = &fbl->soft_buffers[(int) mask];

    file_position_type *filepointer;
    sax_type *saxpointer;

    int current_buffer_number;
    // Check if this buffer is initialized


    if (!current_buffer->initialized)
    {
        pthread_mutex_lock(lock_firstnode);
        if (!current_buffer->initialized)
        {

            current_buffer->max_buffer_size = malloc(sizeof(int)*total_workernumber);
            current_buffer->buffer_size = malloc(sizeof(int)*total_workernumber);
            current_buffer->sax_records=malloc(sizeof(sax_type *)*total_workernumber);
            current_buffer->pos_records=malloc(sizeof(file_position_type *)*total_workernumber);
            for (int i = 0; i < total_workernumber; i++)
            {
                current_buffer->max_buffer_size[i]=0;
                current_buffer->buffer_size[i]=0;
                current_buffer->pos_records[i]=NULL;
                current_buffer->sax_records[i]=NULL;
            }
            current_buffer->node = isax_root_node_init(mask, index->settings->initial_leaf_buffer_size);
            current_buffer->initialized = 1;
            if(index->first_node == NULL)
            {
                index->first_node = current_buffer->node;
                pthread_mutex_unlock(lock_firstnode);
                current_buffer->node->next = NULL;
                current_buffer->node->previous = NULL;

            }
            else
            {
                isax_node * prev_first = index->first_node;
                index->first_node = current_buffer->node;
                index->first_node->next = prev_first;
                prev_first->previous = current_buffer->node;
                pthread_mutex_unlock(lock_firstnode);
            }
            __sync_fetch_and_add(&(index->root_nodes),1);
        }
        else
        {
           pthread_mutex_unlock(lock_firstnode);
        }
    }

    // Check if this buffer is not full!
    if (current_buffer->buffer_size[workernumber] >= current_buffer->max_buffer_size[workernumber]) {
        if(current_buffer->max_buffer_size[workernumber] == 0) {
            current_buffer->max_buffer_size[workernumber] = fbl->initial_buffer_size;
            current_buffer->sax_records[workernumber] = malloc(index->settings->sax_byte_size *
                                                 current_buffer->max_buffer_size[workernumber]);
            current_buffer->pos_records[workernumber] = malloc(index->settings->position_byte_size*
                                                 current_buffer->max_buffer_size[workernumber]);
        }
        else {
            current_buffer->max_buffer_size[workernumber] *= BUFFER_REALLOCATION_RATE;

            current_buffer->sax_records[workernumber] = realloc(current_buffer->sax_records[workernumber],
                                           index->settings->sax_byte_size *
                                           current_buffer->max_buffer_size[workernumber]);
            current_buffer->pos_records[workernumber] = realloc(current_buffer->pos_records[workernumber],
                                           index->settings->position_byte_size *
                                           current_buffer->max_buffer_size[workernumber]);

        }
    }

    if (current_buffer->sax_records[workernumber] == NULL || current_buffer->pos_records[workernumber] == NULL) {
        fprintf(stderr, "error: Could not allocate memory in FBL.");
        return OUT_OF_MEMORY_FAILURE;
    }
    // Copy data to hard buffer and make current buffer point to the hard one
    current_buffer_number=current_buffer->buffer_size[workernumber];
    filepointer=(file_position_type *)current_buffer->pos_records[workernumber];
    saxpointer=(sax_type *)current_buffer->sax_records[workernumber];
    memcpy((void *) (&saxpointer[current_buffer_number*index->settings->n_segments]), (void *) sax, index->settings->sax_byte_size);
    memcpy((void *) (&filepointer[current_buffer_number]), (void *) pos, index->settings->position_byte_size);


    #ifdef DEBUG
    printf("*** Added to node ***\n\n");
    #ifdef TOY
    sax_print(sax, index->settings->n_segments,
              index->settings->sax_bit_cardinality);
    #endif
    #endif
    (current_buffer->buffer_size[workernumber])++;

    return current_buffer->node;
}
static enum response indexconstruction(first_buffer_layer *fbl, isax_index *index,
                                pthread_mutex_t *lock_index, int calculate_thread)
{
    #ifdef DEBUG
    printf("*** FLUSHING ***\n\n");
    #else
    #if VERBOSE_LEVEL == 2
    printf("\n");
    fflush(stdout);
    int i=1;
    fprintf(stdout, "\x1b[31mFlushing: \x1b[36m00.00%%\x1b[0m");
    #endif
    #if VERBOSE_LEVEL == 1
    printf("\n");
    fflush(stdout);
    int i=1;
    fprintf(stdout, "Flushing...\n");
    #endif
    #endif
    trans_fbl_input input_data;
    pthread_t threadid[calculate_thread];
    pthread_mutex_t lock_fbl_conter=PTHREAD_MUTEX_INITIALIZER;
    input_data.index=index;
    input_data.lock_index=lock_index;
    input_data.lock_fbl_conter=&lock_fbl_conter;
    input_data.conternumber=0;
    input_data.stop_number=fbl->number_of_buffers;
    COUNT_QUEUE_TIME_START
    for (int k = 0; k < calculate_thread; k++)
    {
        pthread_create(&(threadid[k]),NULL,indexconstructionworker,(void*)&(input_data));
    }
    for (int k = 0; k < calculate_thread; k++)
    {
        pthread_join(threadid[k],NULL);
    }

    fbl->current_record_index = 0;
    fbl->current_record = fbl->hard_buffer;
    #ifndef DEBUG
    #if VERBOSE_LEVEL == 2
    printf("\n");
    #endif
    #endif

    return SUCCESS;
}
static enum response indexflush(first_buffer_layer *fbl, isax_index *index,
                         pthread_mutex_t *lock_index, int calculate_thread)
{
    #ifdef DEBUG
    printf("*** FLUSHING ***\n\n");
    #else
    #if VERBOSE_LEVEL == 2
    printf("\n");
    fflush(stdout);
    int i=1;
    fprintf(stdout, "\x1b[31mFlushing: \x1b[36m00.00%%\x1b[0m");
    #endif
    #if VERBOSE_LEVEL == 1
    printf("\n");
    fflush(stdout);
    int i=1;
    fprintf(stdout, "Flushing...\n");
    #endif
    #endif
    trans_fbl_input input_data;
    pthread_t threadid[calculate_thread];
    pthread_mutex_t lock_fbl_conter=PTHREAD_MUTEX_INITIALIZER;
    input_data.index=index;
    input_data.lock_index=lock_index;
    input_data.lock_fbl_conter=&lock_fbl_conter;
    input_data.conternumber=0;
    input_data.stop_number=fbl->number_of_buffers;
    COUNT_QUEUE_TIME_START
    for (int k = 0; k < calculate_thread; k++)
    {
        pthread_create(&(threadid[k]),NULL,indexflushworker,(void*)&(input_data));
    }
    for (int k = 0; k < calculate_thread; k++)
    {
        pthread_join(threadid[k],NULL);
    }


    fbl->current_record_index = 0;
    fbl->current_record = fbl->hard_buffer;
    #ifndef DEBUG
    #if VERBOSE_LEVEL == 2
    printf("\n");
    #endif
    #endif

    return SUCCESS;
}
static void *indexconstructionworker(void *input)
{
    isax_index *index=((trans_fbl_input*)input)->index;
    pthread_mutex_t *lock_index=((trans_fbl_input*)input)->lock_index;
    int j,c=1;
    isax_node_record *r = calloc(1, sizeof(isax_node_record));
    while(1)
    {
        pthread_mutex_lock(((trans_fbl_input*)input)->lock_fbl_conter);
        j=((trans_fbl_input*)input)->conternumber;
        ((trans_fbl_input*)input)->conternumber++;
        pthread_mutex_unlock(((trans_fbl_input*)input)->lock_fbl_conter);
        if(j>=((trans_fbl_input*)input)->stop_number)
        {
            break;
        }
        fbl_soft_buffer *current_fbl_node = &index->fbl->soft_buffers[j];
        if (!current_fbl_node->initialized) {
            continue;
        }

        #ifndef DEBUG
        #if VERBOSE_LEVEL == 2
        fprintf(stdout,"\r\x1b[31mFlushing: \x1b[36m%2.2lf%%\x1b[0m", ((float)c/(float)index->root_nodes)*100);
        c++;
        fflush(stdout);
        #endif
        #endif
        int i;
        if (current_fbl_node->buffer_size > 0)
        {
            // For all records in this buffer
            for (i=0; i<current_fbl_node->buffer_size; i++)
            {
                r->sax = (sax_type *) current_fbl_node->sax_records[i];
                r->position = (file_position_type *) current_fbl_node->pos_records[i];
                r->insertion_mode = NO_TMP | PARTIAL;
                // Add record to index
                add_record_to_node(index, current_fbl_node->node, r, 1, 1);

            }

            // flush index node
            flush_subtree_leaf_buffers_m(index, current_fbl_node->node, lock_index);
            // clear FBL records moved in LBL buffers
            free(current_fbl_node->sax_records);
            free(current_fbl_node->pos_records);
            // clear records read from files (free only prev sax buffers)

            isax_index_clear_node_buffers(index, current_fbl_node->node,
                                          INCLUDE_CHILDREN,
                                          TMP_AND_TS_CLEAN);

            index->allocated_memory = 0;
            // Set to 0 in order to re-allocate original space for buffers.
            current_fbl_node->buffer_size = 0;
            current_fbl_node->max_buffer_size = 0;
        }
    }
    free(r);
    return NULL;
}
static void *indexflushworker(void *input)
{
    isax_index *index=((trans_fbl_input*)input)->index;
    pthread_mutex_t *lock_index=((trans_fbl_input*)input)->lock_index;
    int j,c=1;
    while(1)
    {
        pthread_mutex_lock(((trans_fbl_input*)input)->lock_fbl_conter);
        j=((trans_fbl_input*)input)->conternumber;
        ((trans_fbl_input*)input)->conternumber++;
        pthread_mutex_unlock(((trans_fbl_input*)input)->lock_fbl_conter);
        if(j>=((trans_fbl_input*)input)->stop_number)
        {
            break;
        }
        fbl_soft_buffer *current_fbl_node = &index->fbl->soft_buffers[j];
        if (!current_fbl_node->initialized) {
            continue;
        }

        #ifndef DEBUG
        #if VERBOSE_LEVEL == 2
        fprintf(stdout,"\r\x1b[31mFlushing: \x1b[36m%2.2lf%%\x1b[0m", ((float)c/(float)index->root_nodes)*100);
        c++;
        fflush(stdout);
        #endif
        #endif
        if (current_fbl_node->buffer_size > 0)
        {
            // For all records in this buffer

            flush_subtree_leaf_buffers_m(index, current_fbl_node->node, lock_index);

            free(current_fbl_node->sax_records);
            free(current_fbl_node->pos_records);
            // clear records read from files (free only prev sax buffers)

            isax_index_clear_node_buffers(index, current_fbl_node->node,
                                          INCLUDE_CHILDREN,
                                          TMP_AND_TS_CLEAN);

            index->allocated_memory = 0;
            // Set to 0 in order to re-allocate original space for buffers.
            current_fbl_node->buffer_size = 0;
            current_fbl_node->max_buffer_size = 0;
        }
    }
    return NULL;
}

static root_mask_type isax_fbl_index_insert_m(isax_index *index, sax_type *sax,
                                              file_position_type *pos,
                                              pthread_mutex_t *lock_cbl,
                                              pthread_mutex_t *lock_firstnode)
{
    root_mask_type first_bit_mask = 0x00;

    CREATE_MASK(first_bit_mask, index, sax);

    pthread_mutex_lock(&(lock_cbl[first_bit_mask%LOCK_SIZE]));
    insert_to_fbl_m(index->fbl, sax, pos, first_bit_mask, index, lock_firstnode);
    pthread_mutex_unlock(&(lock_cbl[first_bit_mask%LOCK_SIZE]));

    return first_bit_mask;
}
