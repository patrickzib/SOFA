//  
//  sfa_trie.c
//  sfa  C version
//
//  Created by Karima Echihabi on 18/11/2017
//  Copyright 2017 Paris Descartes University. All rights reserved.
//  
//  Based on isax code written by Zoumpatianos on 3/12/12.
//  Copyright 2012 University of Trento. All rights reserved.
//

#include "../config.h"
#include "../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#ifdef VALUES
#include <values.h>
#endif

#include <sys/stat.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "ads/isax_index.h"
#include "ads/dft.h"

#include "ads/sfa.h"

/*
  This functions allocates a two dimensional array
  of num_words rows and (num_symbols-1) columns
  The array will contain the discretization
  intervals
*/
enum response sfa_bins_init(isax_index *index)
{
    int timeseries_size = index->settings->timeseries_size;
    int num_symbols = index->settings->sax_alphabet_cardinality;
    int paa_segments = index->settings->paa_segments;

    index->bins = NULL;  
    index->bins = (ts_type**) calloc(paa_segments,sizeof(ts_type*));
    if(index == NULL) {
        fprintf(stderr,"Error in sfa.c: Could not allocate memory for bins structure.\n");
        return FAILURE;
    }

    //allocate num_symbols-1 memory slots for each word
    for (int i = 0; i < paa_segments; ++i)
    {
        index->bins[i] = calloc(num_symbols-1,sizeof(ts_type));
        if(index == NULL) {
            fprintf(stderr,"Error in sfa.c: Could not allocate memory for bins structure.\n");
            return FAILURE;
        }
        for (int j = 0; j < num_symbols-1; ++j)
        {
            index->bins[i][j] = FLT_MAX;
        }

    }

    fprintf(stderr,"Initialized bins[%d][%d] \n", paa_segments, num_symbols-1 );      
    return SUCCESS;
}

void sfa_free_bins(isax_index *index)
{
    for (int i = 0; i < index->settings->paa_segments; ++i)
    {
        free(index->bins[i]);
    }
    free(index->bins);
}

void sfa_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread)
{
    int paa_segments = index->settings->paa_segments;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;
    
    fprintf(stderr, ">>> Binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START

    ts_type ** dft_mem_array = (ts_type **) calloc(paa_segments, sizeof(ts_type*));
    for (int k = 0; k < paa_segments; ++k)
    {
        dft_mem_array[k] = (ts_type *) calloc(sample_size,sizeof(ts_type));
    }

    //first build the bins out of a sample of the data
    //read whole sample in memory
    pthread_t threadid[maxquerythread];
    bins_data_inmemory *input_data=malloc(sizeof(bins_data_inmemory)*(maxquerythread));

    ts_type * ts;
    fftwf_complex *ts_out;
    fftwf_plan plan_forward;
    ts_type * transform;
    
    for (int i = 0; i < maxquerythread; i++)
    {
        ts = fftwf_malloc ( sizeof ( ts_type ) * ts_length);
        ts_out = (fftwf_complex *)fftwf_malloc ( sizeof ( fftwf_complex ) * (ts_length/2+1) );
        plan_forward = fftwf_plan_dft_r2c_1d (ts_length, ts, ts_out, FFTW_ESTIMATE);
        transform = fftwf_malloc ( sizeof ( ts_type ) * ts_length);

        input_data[i].index=index;
        input_data[i].dft_mem_array=dft_mem_array;
        input_data[i].filename=ifilename;
        input_data[i].workernumber=i;
        input_data[i].records = sample_size/maxquerythread;
        input_data[i].records_offset = sample_size/maxquerythread;

        //first-n-sampling
        if(index->settings->sample_type == 1)
        {
            input_data[i].start_number=i*(sample_size/maxquerythread);
            input_data[i].stop_number=(i+1)*(sample_size/maxquerythread);
        }
        //uniform sampling
        else if(index->settings->sample_type == 2)
        {
            input_data[i].start_number=i*(ts_num/maxquerythread);
            input_data[i].stop_number=(i+1)*(ts_num/maxquerythread);
        }
        //random sampling
        else if(index->settings->sample_type == 3)
        {
            input_data[i].start_number=0;
            input_data[i].stop_number=ts_num;
        }
        
        input_data[i].ts = ts;
        input_data[i].ts_out = ts_out;
        input_data[i].plan_forward = plan_forward;
        input_data[i].transform = transform;
    }

    //reset values for last worker to keep lost segments at the end
    input_data[maxquerythread-1].records = sample_size - (maxquerythread-1)*(sample_size/maxquerythread);

    if(index->settings->sample_type == 1)
    {
        input_data[maxquerythread-1].stop_number=sample_size;
    }
    else if(index->settings->sample_type == 2)
    {
        input_data[maxquerythread-1].stop_number=ts_num;
    }

    for (int i = 0; i < maxquerythread; i++)
    {
        pthread_create(&(threadid[i]),NULL,set_bins_worker_dft,(void*)&(input_data[i]));
    }

    //wait for the finish of other threads
    for (int i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
    }
  
    for (int i = 0; i < maxquerythread; i++)
    {   
        fftwf_destroy_plan (input_data[i].plan_forward);
        fftwf_free (input_data[i].ts);
        fftwf_free (input_data[i].ts_out);
        fftwf_free (input_data[i].transform);

        input_data[i].start_number=i*(paa_segments/maxquerythread);
        input_data[i].stop_number=(i+1)*(paa_segments/maxquerythread);
    }

    input_data[maxquerythread-1].start_number=(maxquerythread-1)*(paa_segments/maxquerythread);
    input_data[maxquerythread-1].stop_number=paa_segments;

    for (int i = 0; i < maxquerythread; i++)
    {
        pthread_create(&(threadid[i]),NULL,order_divide_worker,(void*)&(input_data[i]));
    }
    //wait for the finish of other threads
    for (int i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
    }

    free(input_data);
    free_dft_memory(index, dft_mem_array);

    COUNT_BINNING_TIME_END

    sfa_print_bins(index);
	fprintf(stderr, ">>> Finished binning\n");
}

void* set_bins_worker_dft(void *transferdata)
{
    struct bins_data_inmemory* bins_data = (bins_data_inmemory*)transferdata;

    ts_type ** dft_mem_array = bins_data->dft_mem_array;

    isax_index *index= ((bins_data_inmemory*)transferdata)->index;
    unsigned long start_number=bins_data->start_number;
    unsigned long stop_number=bins_data->stop_number;

    unsigned long ts_length = index->settings->timeseries_size;
    int paa_segments =  index->settings->paa_segments;
    
    unsigned long start_index = start_number*ts_length*sizeof(ts_type);

    FILE *ifile;
    ifile = fopen (bins_data->filename,"rb");
    fseek (ifile , start_index, SEEK_SET);

    unsigned long skip_elements;
    long records = bins_data->records;

    //set number of elements to skip for uniform sampling
    if(index->settings->sample_type == 2)
    {
        skip_elements = (((stop_number-start_number)/records)-1)*ts_length*sizeof(ts_type);
    }

    ts_type *ts = bins_data->ts;
    fftwf_complex *ts_out=bins_data->ts_out;
    fftwf_plan plan_forward=bins_data->plan_forward;   
  	ts_type* transform = bins_data->transform;

  	ts_type * ts_orig = (ts_type *) calloc(index->settings->timeseries_size, sizeof(ts_type));

    for(int i=0; i<records; ++i)
    {
        //choose random position for random sampling
        if(index->settings->sample_type == 3)
        {
            long int position = (rand()%(stop_number-start_number+1)+start_number);
            fseek(ifile, (position*ts_length*sizeof(ts_type)), SEEK_SET);
        }

        fread(ts_orig, sizeof(ts_type), ts_length, ifile);
            
        for (int j =0; j< ts_length; ++j)
        {
            ts[j] = ts_orig[j];
        }

        fft_from_ts(index, ts, ts_out, transform, plan_forward);
        
        for (int j = 0; j < paa_segments; ++j)      
        {
            ts_type value =(ts_type) roundf(transform[j]*100.0)/100.0;
            dft_mem_array[j][i + (bins_data->workernumber * bins_data->records_offset)]=value;   
        }

        //skip elements for uniform sampling
        if(index->settings->sample_type == 2)
        {
            fseek(ifile, skip_elements, SEEK_CUR);
        }
    }

    free(ts_orig);
    fclose(ifile);
}

void* order_divide_worker(void *transferdata)
{
    struct bins_data_inmemory* bins_data = (bins_data_inmemory*)transferdata;

    ts_type ** dft_mem_array = bins_data->dft_mem_array;

    isax_index *index= ((bins_data_inmemory*)transferdata)->index;
    unsigned long start_number=bins_data->start_number;
    unsigned long stop_number=bins_data->stop_number;

    unsigned int sample_size = index->settings->sample_size;
    int paa_segments =  index->settings->paa_segments;
    ts_type * cur_coeff_line;

    for (int j = start_number; j < stop_number; ++j)
    {
        cur_coeff_line = (ts_type *) dft_mem_array[j];   
        qsort(cur_coeff_line, sample_size, sizeof(ts_type), &compare_ts_type);
    }
    //equi-depth splitting
    if(index->settings->histogram_type==1)
    {
        ts_type depth = (ts_type) sample_size / index->settings->sax_alphabet_cardinality;

        for (int i = start_number; i < stop_number; ++i)
        {
            cur_coeff_line = dft_mem_array[i];

            int pos = 0;
            int count = 0;
            for (int j=0; j < sample_size; ++j)
            {
                ++count;
                unsigned int ceiling = ceil (depth * (pos+1));
                if (count > ceiling && (pos == 0 || index->bins[i][pos-1] != cur_coeff_line[j]))
                {
                    index->bins[i][pos] = cur_coeff_line[j];
                    pos++;
                }
            }
        }
    }
    //equi-width splitting
    else if(index->settings->histogram_type==2)
    {
        int num_symbols = index->settings->sax_alphabet_cardinality;

        for (int i = start_number; i < stop_number; ++i)
        {
            cur_coeff_line = dft_mem_array[i];
           
            ts_type first = cur_coeff_line[0];
            ts_type last = cur_coeff_line[sample_size-1];

            ts_type interval_width = (last-first) / (ts_type) num_symbols;
            for (int j=0; j < num_symbols-1; ++j)
            {
                index->bins[i][j] = interval_width*(j+1)+first;
            }
        }
    }
}
//deprecated
void sfa_fill_order_line(isax_index *index, ts_type **dft_mem_array)
{
    unsigned int sample_size = index->settings->sample_size;

    int transforms_size = index->settings->paa_segments;

    fprintf(stderr,"Sample size %u\n", sample_size);     

    ts_type * cur_coeff_line;

    for (int j = 0; j < transforms_size; ++j)
    {
        cur_coeff_line = (ts_type *) dft_mem_array[j];      
        qsort(cur_coeff_line, sample_size, sizeof(ts_type), &compare_ts_type);
    }
}

//deprecated
void sfa_divide_equi_depth_hist(isax_index *index, ts_type **dft_mem_array)
{
    unsigned int sample_size = index->settings->sample_size;
    ts_type * cur_coeff_line;

    ts_type depth = (ts_type) sample_size / index->settings->sax_alphabet_cardinality;
    fprintf(stderr, "Using Equi-depth histograms\n");

    for (int i = 0; i < index->settings->paa_segments; ++i)
    {
        cur_coeff_line = dft_mem_array[i];

	    int pos = 0;
	    int count = 0;
        for (int j=0; j < sample_size; ++j)
	    {
	        ++count;
	        unsigned int ceiling = ceil (depth * (pos+1));
	        if (count > ceiling && (pos == 0 || index->bins[i][pos-1] != cur_coeff_line[j]))
	        {
	            index->bins[i][pos] = cur_coeff_line[j];
	            pos++;
 	        }
	    }
    }
}
//deprecated
void sfa_divide_equi_width_hist(isax_index *index, ts_type **dft_mem_array)
{
    ts_type * cur_coeff_line;
    unsigned int sample_size = index->settings->sample_size;
    int num_symbols = index->settings->sax_alphabet_cardinality;

    fprintf(stderr, "Using Equi-width histograms\n");
  
    for (int i = 0; i < index->settings->paa_segments; ++i)
    {
        cur_coeff_line = dft_mem_array[i];
	   
        ts_type first = cur_coeff_line[0];
        ts_type last = cur_coeff_line[sample_size-1];

        ts_type interval_width = (last-first) / (ts_type) num_symbols;
        for (int j=0; j < num_symbols-1; ++j)
        {
	        index->bins[i][j] = interval_width*(j+1)+first;
        }
    }
}

void sfa_print_bins(isax_index *index)
{
    fprintf(stderr,"Sample size %u\n", index->settings->sample_size);     
	if(index->settings->histogram_type==1)
	{
		fprintf(stderr, "Using Equi-depth histograms\n");
	}
	else if(index->settings->histogram_type==2)
	{
		fprintf(stderr, "Using Equi-width histograms\n");
	}
    int paa_segments = index->settings->paa_segments;
    fprintf(stderr,"[\n");
    for (int i = 0; i < paa_segments; ++i)
    {
        fprintf(stderr,"-Inf\t");
        for (int j=0; j < index->settings->sax_alphabet_cardinality-1; ++j)
        {
            ts_type value = roundf(index->bins[i][j]*100.0)/100.0;
            if (value == FLT_MAX)	  
	            fprintf(stderr,", Inf\t\n");
	        else
	            fprintf(stderr,",\t%g\t",value);
        }
        fprintf(stderr,";\n");
    }
    fprintf(stderr,"]\n");
}

void free_dft_memory(isax_index *index, ts_type **dft_mem_array)
{
    int paa_segments = index->settings->paa_segments;
    for (int k = 0; k < paa_segments; ++k)
    {
        free(dft_mem_array[k]);
    }        
    free(dft_mem_array);
}

int compare_ts_type (const void * a, const void * b)
{
    
    ts_type ts_a = *((ts_type*) a);
    ts_type ts_b = *((ts_type*) b);
    
    if (ts_a < ts_b )
        return -1;
    else return ts_a>ts_b;
}

ts_type sfa_fft_min_dist (isax_index *index, unsigned char c1_value, unsigned char c2_value, ts_type real_c2, unsigned int dim)
{
    if (c1_value == c2_value) {
        return 0;
    }

    if (c1_value > c2_value)
    {
        return  (index->bins[dim][((int)c1_value)-1] - real_c2);
    }
  
    if (c1_value < c2_value)
    {
        return  (real_c2 - index->bins[dim][(int)c1_value]);
    }
}

ts_type minidist_fft_to_isax_complete(isax_index *index, ts_type *query_fft, sax_type *sax, sax_type *sax_cardinalities, float min_val) 
{ 
    sax_type * query_sax = calloc(index->settings->paa_segments,sizeof(sax_type));

    sfa_from_fft(index, (ts_type *) query_fft, query_sax);
    ts_type distance = 0.0;
    ts_type value;

    unsigned int i=0;

    if (!index->settings->is_norm)
    {
        sax_type sax_promoted = sax[0];

        //promote sax , if necessary
        if((unsigned int) sax_cardinalities[0] < index->settings->sax_bit_cardinality)
        {
           sax_promoted = promote_sax_cardinality(index, query_fft[0], query_sax[0], sax[0], sax_cardinalities[0]);
        }

        //distance = sfa_fft_min_dist(index, sax_promoted, sax[0], query_fft[0],0);
        if (sax_promoted == query_sax[0]) 
        {
            distance = 0.0;
        }
        else if (sax_promoted > query_sax[0])
        {
            distance = (index->bins[0][((int)sax_promoted)-1] - query_fft[0]);
        }
        else
        {
            distance = (query_fft[0] - index->bins[0][(int)sax_promoted]);
        }

        distance *= distance;
        i += 2;
    }

    for(;i<index->settings->paa_segments; ++i)
    {

        sax_type sax_promoted = sax[i];

        //promote sax , if necessary
        if((unsigned int) sax_cardinalities[i] < index->settings->sax_bit_cardinality)
        {
           sax_promoted = promote_sax_cardinality(index, query_fft[i], query_sax[i], sax[i], sax_cardinalities[i]);
        }

        //distance = sfa_fft_min_dist(index, sax_promoted, sax[i], query_fft[i],i);
        if (sax_promoted == query_sax[i]) 
        {
            value = 0.0;
        }
        else if (sax_promoted > query_sax[i])
        {
            value = (index->bins[i][((int)sax_promoted)-1] - query_fft[i]);
        }
        else
        {
            value = (query_fft[i] - index->bins[i][(int)sax_promoted]);
        }

        distance += 2*value*value;

        if(distance>min_val)
        {
            free(query_sax);
            return distance;
        }
    }

    free(query_sax);
    return distance;
}

sax_type promote_sax_cardinality(isax_index *index, ts_type current_query_fft, sax_type current_query_sax, sax_type current_value_sax, sax_type current_sax_cardinality)
{
    if((unsigned int) current_sax_cardinality < index->settings->sax_bit_cardinality)
    {
        current_value_sax = current_value_sax << (index->settings->sax_bit_cardinality - current_sax_cardinality);

        sax_type mask_bits = (0xff << (index->settings->sax_bit_cardinality - current_sax_cardinality));
        sax_type current_query_sax_masked = current_query_sax & mask_bits;

        //lex. smaller -> set all  other values 1
        if(current_value_sax < current_query_sax_masked)
        {
            mask_bits = (0xff >> current_sax_cardinality);
            return current_value_sax | mask_bits;
        }
        //lex. larger -> set all other values 0
        else if(current_value_sax > current_query_sax_masked)
        {
            return current_value_sax;
        }
        //lex. same -> set all other values like query
        else if(current_value_sax == current_query_sax_masked)
        {
            return current_query_sax;
        }
    }
    else if((unsigned int) current_sax_cardinality == index->settings->sax_bit_cardinality)
    {
        return current_value_sax;
    }
    else if((unsigned int) current_sax_cardinality > index->settings->sax_bit_cardinality)
    {
        fprintf(stderr, "ERROR: current cardinality %u is bigger than maximal cardinality %u\n",(unsigned int) current_sax_cardinality, index->settings->sax_bit_cardinality); 
    }
}

ts_type minidist_fft_to_isax(isax_index *index, ts_type *query_fft, sax_type *sax, float min_val) 
{
    ts_type distance = 0.0;
    ts_type value;
    unsigned int i=0;

    sax_type * query_sax = calloc(index->settings->paa_segments,sizeof(sax_type));

    sfa_from_fft(index, (ts_type *) query_fft, query_sax);

    if (!index->settings->is_norm)
    {
        //distance = sfa_fft_min_dist(index, sax[0], query_sax[0], query_fft[0],0);
        if (sax[0] == query_sax[0]) 
        {
            distance = 0.0;
        }
        else if (sax[0] > query_sax[0])
        {
            distance = (index->bins[0][((int)sax[0])-1] - query_fft[0]);
        }
        else
        {
            distance = (query_fft[0] - index->bins[0][(int)sax[0]]);
        }

        distance *= distance;

        i += 2;
    }

    for (; i < index->settings->paa_segments; i++) {

        //ts_type value = sfa_fft_min_dist(index, sax[i], query_sax[i], query_fft[i],i);
        if (sax[i] == query_sax[i]) 
        {
            value = 0.0;
        }
        else if (sax[i] > query_sax[i])
        {
            value = (index->bins[i][((int)sax[i])-1] - query_fft[i]);
        }
        else
        {
            value = (query_fft[i] - index->bins[i][(int)sax[i]]);
        }

        distance += 2*value*value;

        if(distance>min_val)
        {
            free(query_sax);
            return distance;
        }
    }
    free(query_sax);
    return distance;
}

void sfa_printbin(unsigned long long n, int size)
{
    char *b = malloc(sizeof(char) * (size+1));
    int i;
    
    for (i=0; i<size; i++) {
        b[i] = '0';
    }
    
    for (i=0; i<size; i++, n=n/2)
        if (n%2) b[size-1-i] = '1';
    
    b[size] = '\0';
    printf("%s\n", b);
    free(b);
}

void sfa_print(sax_type *sax, int segments, int cardinality) 
{
    int i;
    for (i=0; i < segments; i++) {
        printf("%d:\t\n", i);
        sfa_printbin(sax[i], cardinality);
    }
    printf("\n");
}

void fft_print(ts_type *fft, int segments)
{
    int i;
    for (i=0; i < segments; i++) {
        printf("%d:\t%.3f\n", i, fft[i]);
    }
    printf("\n");
}