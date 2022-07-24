//  
//  sfa.c
//  sfa  C version for MESSI
//
//  Based on sfa_trie code by Karima Echihabi on 18/11/2017
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
#include "immintrin.h"
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
    index->binsv = (ts_type*) calloc(paa_segments*(num_symbols-1),sizeof(ts_type));
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
    for (int j = 0; j < paa_segments*(num_symbols-1); ++j)
    {
        index->binsv[j] = FLT_MAX;
    }
    fprintf(stderr,"Initialized bins[%d][%d] \n", paa_segments, num_symbols-1 );

    if(index->settings->coeff_number != 0)
    {
        index->coefficients = calloc(paa_segments/2, sizeof(int));
    } 

    return SUCCESS;
}

/*
  This functions frees the allocated bins-array
*/
void sfa_free_bins(isax_index *index)
{
    for (int i = 0; i < index->settings->paa_segments; ++i)
    {
        free(index->bins[i]);
    }
    free(index->bins);
}

/*
  In this function, the intervals are caluclated (multiple coeff. binning) and saved to bins
*/
void sfa_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread)
{
    int paa_segments = index->settings->paa_segments;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;
    
    fprintf(stderr, ">>> Binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START

    //alocate dft_mem_array for to savetime series from sampling
    ts_type ** dft_mem_array = (ts_type **) calloc(paa_segments, sizeof(ts_type*));
    for (int k = 0; k < paa_segments; ++k)
    {
        dft_mem_array[k] = (ts_type *) calloc(sample_size,sizeof(ts_type));
    }

    //build the bins out of a sample of the data
    //read whole sample in memory
    pthread_t threadid[maxquerythread];
    bins_data_inmemory *input_data=malloc(sizeof(bins_data_inmemory)*(maxquerythread));

    ts_type * ts;
    fftwf_complex *ts_out;
    fftwf_plan plan_forward;
    ts_type * transform;
    
    for (int i = 0; i < maxquerythread; i++)
    {
        //create FFTW-objects for the threads
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

    //initiate worker threads for sampling values and calculating FFTs
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

    //initiate worker threads for splitting coefficients into intervals
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

/*
  In this function, the intervals are caluclated (multiple coeff. binning).
  The coefficients with the highest variance in the intervals are chosen and these values are saved to bins
*/
void sfa_set_bins_coeff(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread)
{
    int paa_segments = index->settings->paa_segments;
    int coeff_number = index->settings->coeff_number;
    int ts_length = index->settings->timeseries_size;
    unsigned int sample_size = index->settings->sample_size;
    
    fprintf(stderr, ">>> Binning: %s\n", ifilename);
    COUNT_BINNING_TIME_START

    ts_type ** dft_mem_array = (ts_type **) calloc(coeff_number, sizeof(ts_type*));
    for (int k = 0; k < coeff_number; ++k)
    {
        dft_mem_array[k] = (ts_type *) calloc(sample_size,sizeof(ts_type));
    }

    //build the bins out of a sample of the data
    //read whole sample in memory
    pthread_t threadid[maxquerythread];
    bins_data_inmemory *input_data=malloc(sizeof(bins_data_inmemory)*(maxquerythread));

    ts_type * ts;
    fftwf_complex *ts_out;
    fftwf_plan plan_forward;
    ts_type * transform;
    
    for (int i = 0; i < maxquerythread; i++)
    {
        //create FFTW-objects for the threads
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

    //initiate worker threads for sampling values and calculating FFTs
    for (int i = 0; i < maxquerythread; i++)
    {
        pthread_create(&(threadid[i]),NULL,set_bins_worker_dft_coeff,(void*)&(input_data[i]));
    }

    //wait for the finish of other threads
    for (int i = 0; i < maxquerythread; i++)
    {
        pthread_join(threadid[i],NULL);
    }

    //calculate coefficient-wise variance
    ts_type **dft_mem_array_coeff = calculate_variance_coeff(index, dft_mem_array);

    for(int i=0; i<coeff_number; ++i)
    {
        free(dft_mem_array[i]);
    }
    free(dft_mem_array);

    for (int i = 0; i < maxquerythread; i++)
    {   
        fftwf_destroy_plan (input_data[i].plan_forward);
        fftwf_free (input_data[i].ts);
        fftwf_free (input_data[i].ts_out);
        fftwf_free (input_data[i].transform);

        input_data[i].dft_mem_array=dft_mem_array_coeff;

        input_data[i].start_number=i*(paa_segments/maxquerythread);
        input_data[i].stop_number=(i+1)*(paa_segments/maxquerythread);
    }

    input_data[maxquerythread-1].start_number=(maxquerythread-1)*(paa_segments/maxquerythread);
    input_data[maxquerythread-1].stop_number=paa_segments;

    //initiate worker threads for splitting coefficients into intervals
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
    free_dft_memory(index, dft_mem_array_coeff);

    COUNT_BINNING_TIME_END

    sfa_print_bins(index);
    fprintf(stderr, ">>> Finished binning\n");
}

/*
  This function calculates the variance for each coefficient in dft_mem_array.
  It returns a trimmed dft_mem_array with only the highest-variance coeff.
*/
ts_type** calculate_variance_coeff(isax_index *index, ts_type ** dft_mem_array)
{
    int coeff_number = index->settings->coeff_number;
    int paa_segments = index->settings->paa_segments;
    unsigned int sample_size = index->settings->sample_size;

    struct variance_coeff_index var_coeff_index[coeff_number/2];
    
    for(int i=0; i<coeff_number/2; ++i)
    {
        double mean_real = 0.0;
        double mean_imag = 0.0;
        double var_real = 0.0;
        double var_imag = 0.0;

        for(int j=0; j<sample_size; ++j)
        {
            mean_real += dft_mem_array[i*2][j];
            mean_imag += dft_mem_array[i*2+1][j];
        }
        mean_real = mean_real / (double) sample_size;
        mean_imag = mean_imag / (double) sample_size;

        for(int j=0; j<sample_size; ++j)
        {
            var_real += (dft_mem_array[i*2][j] - mean_real)*(dft_mem_array[i*2][j] - mean_real);
            var_imag += (dft_mem_array[i*2+1][j] - mean_imag)*(dft_mem_array[i*2+1][j] - mean_imag);
        }
        var_real = var_real / (double) sample_size;
        var_imag = var_imag / (double) sample_size;

        double total_var = var_real + var_imag;

        var_coeff_index[i].variance = total_var;
        var_coeff_index[i].coeff_index = i;
    }

    for(int i=0; i<coeff_number/2; ++i)
    {
        fprintf(stderr,"variance %.3f\tposition %d\n", var_coeff_index[i].variance, var_coeff_index[i].coeff_index);
    }

    qsort(var_coeff_index, coeff_number/2, sizeof(var_coeff_index[0]), compare_var);

    fprintf(stderr,"SORTED:\n");
    for(int i=0; i<coeff_number/2; ++i)
    {
        fprintf(stderr,"variance %.3f\tposition %d\n", var_coeff_index[i].variance, var_coeff_index[i].coeff_index);
    }

    for(int i=0; i<paa_segments/2; ++i)
    {
        index->coefficients[i] = var_coeff_index[i].coeff_index;
    }
    qsort(index->coefficients, paa_segments/2, sizeof(int), compare_int);

    fprintf(stderr,"HIGHEST COEFFS SORTED:\n");
    for(int i=0; i<paa_segments/2; ++i)
    {
        fprintf(stderr,"%d\n",index->coefficients[i]);
    }

    ts_type ** dft_mem_array_coeff = (ts_type **) calloc(paa_segments, sizeof(ts_type*));
    for (int k = 0; k < paa_segments; ++k)
    {
        dft_mem_array_coeff[k] = (ts_type *) calloc(sample_size,sizeof(ts_type));
    }
    
    for(int i=0; i<paa_segments/2; ++i)
    {
        int coeff = index->coefficients[i];

        for(int j=0; j<sample_size; ++j)
        {
            dft_mem_array_coeff[i*2][j] = dft_mem_array[coeff*2][j];
            dft_mem_array_coeff[i*2+1][j] = dft_mem_array[coeff*2+1][j];
        }
    }

    return dft_mem_array_coeff;
}

/*
    Worker method for sampling values, calculating FFT coefficients and saving them to dft_mem_array
*/
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
            //long int position = (rand()%(stop_number-start_number+1)+start_number);
            long int position = start_number + random_at_most(records);
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
/*
    Worker method for sampling values, calculating FFT coefficients (the first coeff_number coefficients) and saving them to dft_mem_array
*/
void* set_bins_worker_dft_coeff(void *transferdata)
{
    struct bins_data_inmemory* bins_data = (bins_data_inmemory*)transferdata;

    ts_type ** dft_mem_array = bins_data->dft_mem_array;

    isax_index *index= ((bins_data_inmemory*)transferdata)->index;
    unsigned long start_number=bins_data->start_number;
    unsigned long stop_number=bins_data->stop_number;

    unsigned long ts_length = index->settings->timeseries_size;
    int coeff_number = index->settings->coeff_number;
    
    unsigned long start_index = start_number*ts_length*sizeof(ts_type);

    FILE *ifile;
    ifile = fopen (bins_data->filename,"rb");
    fseek (ifile , start_index, SEEK_SET);

    unsigned long skip_elements;
    long records = bins_data->records;

    unsigned long position_count = start_number;

    //set number of elements to skip for uniform sampling
    if(index->settings->sample_type == 2)
    {
        skip_elements = (((stop_number-start_number)/records)-1);
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
            long int position = start_number + random_at_most(records);
            fseek(ifile, (position*ts_length*sizeof(ts_type)), SEEK_SET);
        }

        fread(ts_orig, sizeof(ts_type), ts_length, ifile);
            
        for (int j =0; j< ts_length; ++j)
        {
            ts[j] = ts_orig[j];
        }

        fft_from_ts_all_coeff(index, ts, ts_out, transform, plan_forward);
        
        for (int j = 0; j < coeff_number; ++j)      
        {
            ts_type value =(ts_type) roundf(transform[j]*100.0)/100.0;
            dft_mem_array[j][i + (bins_data->workernumber * bins_data->records_offset)]=value;   
        }

        //skip elements for uniform sampling
        if(index->settings->sample_type == 2)
        {

		fseek(ifile, skip_elements*ts_length*sizeof(ts_type), SEEK_CUR);
		position_count += (1+skip_elements);
        	if(position_count>=stop_number)
        	{
        		fprintf(stderr,"pos %lu; stop_number %lu\n",position_count,stop_number);
        	}
        }
    }

    free(ts_orig);
    fclose(ifile);
}

/*
    Worker method for coefficient-wise splitting and saving the interval to bins
*/
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


/*
    Method for printing the discretization intervals to the console
*/
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

//compare-functions for qsort
int compare_ts_type (const void *a, const void *b)
{
    ts_type ts_a = *((ts_type*) a);
    ts_type ts_b = *((ts_type*) b);
    
    if (ts_a < ts_b )
        return -1;
    else return ts_a>ts_b;
}

int compare_var (const void *a, const void *b)
{
    struct variance_coeff_index *a1 = (struct variance_coeff_index *)a;
    struct variance_coeff_index *a2 = (struct variance_coeff_index *)b;
    if ((*a1).variance > (*a2).variance)
        return -1;
    else if ((*a1).variance < (*a2).variance)
        return 1;
    else
        return 0;
}

int compare_int (const void *a, const void *b)
{
    const int *ia = (const int *)a;
    const int *ib = (const int *)b;
    return *ia  - *ib; 
}

/*
    This function calculates a mindist (lower bounding dist.) between a query (FFT coeff.) and a SFA representation
*/
ts_type minidist_fft_to_isax(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf)
{
    int min_val = MINVAL;
    int max_val = MAXVAL;

    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->paa_segments;
   
    ts_type distance = 0.0;
    ts_type value;

    int i=0;

    //special case: for not normalized time series, the first coefficient has to be treated spacially
    //for normalized data, this part is skipped
    if(!index->settings->is_norm && (index->settings->coeff_number==0 || index->coefficients[0]==0))
    {
        sax_type c_c = sax_cardinalities[i];

        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        
        sax_type region_lower = (v <<  (c_m - c_c));
        sax_type region_upper = (~((int)MAXFLOAT << (c_m - c_c)) | region_lower);

        float breakpoint_lower = 0;
        float breakpoint_upper = 0;
        
        
        if (region_lower == 0) {
            breakpoint_lower = min_val;
        }
        else
        {
            breakpoint_lower = index->bins[i][ region_lower - 1];
        }
        if (region_upper == max_cardinality - 1) {
            breakpoint_upper = max_val;
        }
        else
        {
            breakpoint_upper = index->bins[i][ region_upper];
        }
        
        //real part for first coeff. is added without factor 2
        if (breakpoint_lower > fft[i]) {

            distance = breakpoint_lower - fft[i];
            distance *= distance;
        }
        else if(breakpoint_upper < fft[i]) {
            distance = (fft[i] - breakpoint_upper);
            distance *= distance;
        }

        //skip the imaginary part of the first coefficient
        i = 2;
    }

    for (; i<number_of_segments; i++) {

        sax_type c_c = sax_cardinalities[i];

        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        
        sax_type region_lower = (v <<  (c_m - c_c));
        sax_type region_upper = (~((int)MAXFLOAT << (c_m - c_c)) | region_lower);

        float breakpoint_lower = 0;
        float breakpoint_upper = 0;
        
        
        if (region_lower == 0) {
            breakpoint_lower = min_val;
        }
        else
        {
            breakpoint_lower = index->bins[i][ region_lower - 1];
        }
        if (region_upper == max_cardinality - 1) {
            breakpoint_upper = max_val;
        }
        else
        {
            breakpoint_upper = index->bins[i][ region_upper];
        }
        
        if (breakpoint_lower > fft[i]) {

            value = breakpoint_lower - fft[i];
            distance += 2*value*value;
        }
        else if(breakpoint_upper < fft[i]) {
            value = (fft[i] - breakpoint_upper);
            distance += 2*value*value;
        }

        if(distance>bsf)
        {
            return distance;
        }
    }

    return distance;
}

/*
    This function calculates a mindist (lower bounding dist.) between a query (FFT coeff.) and a SFA representation
*/
ts_type minidist_fft_to_isax_raw(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf)
{
    int min_val = MINVAL;
    int max_val = MAXVAL;

    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->paa_segments;
   
    ts_type distance = 0.0;
    ts_type value;
    
    int i=0;

    //special case: for not normalized time series, the first coefficient has to be treated spacially
    //for normalized data, this part is skipped
    if(!index->settings->is_norm && (index->settings->coeff_number==0 || index->coefficients[0]==0))
    {
        sax_type c_c = sax_cardinalities[i];

        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        
        sax_type region_lower = (v <<  (c_m - c_c));
        sax_type region_upper = (~((int)MAXFLOAT << (c_m - c_c)) | region_lower);

        float breakpoint_lower = 0;
        float breakpoint_upper = 0;
        
        
        if (region_lower == 0) {
            breakpoint_lower = min_val;
        }
        else
        {
            breakpoint_lower = index->bins[i][region_lower - 1];
        }
        if (region_upper == max_cardinality - 1) {
            breakpoint_upper = max_val;
        }
        else
        {
            breakpoint_upper = index->bins[i][ region_upper];
        }
        
        //real part for first coeff. is added without factor 2
        if (breakpoint_lower > fft[i]) {

            distance = breakpoint_lower - fft[i];
            distance *= distance;
        }
        else if(breakpoint_upper < fft[i]) {
            distance = (fft[i] - breakpoint_upper);
            distance *= distance;
        }

        //skip the imaginary part of the first coefficient
        i = 2;
    }

    for (; i<number_of_segments; i++) {

        sax_type c_c = sax_cardinalities[i];

        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        
        sax_type region_lower = (v <<  (c_m - c_c));
        sax_type region_upper = (~((int)MAXFLOAT << (c_m - c_c)) | region_lower);

        float breakpoint_lower = 0;
        float breakpoint_upper = 0;
        
        
        if (region_lower == 0) {
            breakpoint_lower = min_val;
        }
        else
        {
            breakpoint_lower = index->bins[i][  region_lower - 1];
        }
        if (region_upper == max_cardinality - 1) {
            breakpoint_upper = max_val;
        }
        else
        {
            breakpoint_upper = index->bins[i][ region_upper];
        }
        
        if (breakpoint_lower > fft[i]) {

            value = breakpoint_lower - fft[i];
            distance += 2*value*value;
        }
        else if(breakpoint_upper < fft[i]) {
            value = (fft[i] - breakpoint_upper);
            distance += 2*value*value;
        }

        if(distance>bsf)
        {
            return distance;
        }
    }

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

/*
    This function calculates random numbers between 0 and max
*/
long random_at_most(long max) {
	unsigned long
	num_bins = (unsigned long) max + 1,
	num_rand = (unsigned long) RAND_MAX + 1,
	bin_size = num_rand / num_bins,
	defect   = num_rand % num_bins;

	long x;
	do{
		x = random();
	}
	while (num_rand - defect <= (unsigned long)x);

	return x/bin_size;
}
ts_type minidist_fft_to_isax_raw_autoSIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf)
{
    int min_val = MINVAL;
    int max_val = MAXVAL;

    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;
    int max_cardinality = index->settings->sax_alphabet_cardinality;
    int number_of_segments = index->settings->paa_segments;
   
    ts_type distance = 0.0;
    ts_type value;
    
    int i=0;

    //special case: for not normalized time series, the first coefficient has to be treated spacially
    //for normalized data, this part is skipped
    if(!index->settings->is_norm && (index->settings->coeff_number==0 || index->coefficients[0]==0))
    {
        sax_type c_c = sax_cardinalities[i];

        sax_type c_m = max_bit_cardinality;
        sax_type v = sax[i];
        
        sax_type region_lower = (v <<  (c_m - c_c));
        sax_type region_upper = (~((int)MAXFLOAT << (c_m - c_c)) | region_lower);

        float breakpoint_lower = 0;
        float breakpoint_upper = 0;
        
        
        if (region_lower == 0) {
            breakpoint_lower = min_val;
        }
        else
        {
            breakpoint_lower = index->bins[i][region_lower - 1];
        }
        if (region_upper == max_cardinality - 1) {
            breakpoint_upper = max_val;
        }
        else
        {
            breakpoint_upper = index->bins[i][ region_upper];
        }
        
        //real part for first coeff. is added without factor 2
        if (breakpoint_lower > fft[i]) {

            distance = breakpoint_lower - fft[i];
            distance *= distance;
        }
        else if(breakpoint_upper < fft[i]) {
            distance = (fft[i] - breakpoint_upper);
            distance *= distance;
        }

        //skip the imaginary part of the first coefficient
        i = 2;
    }
    sax_type c_cv[number_of_segments];
    sax_type vv[number_of_segments];
    sax_type c_m = max_bit_cardinality;
   for (; i<number_of_segments; i++) 
    {
          c_cv[i]  = sax_cardinalities[i];
          vv[i]= sax[i];
    }
    sax_type region_lowerv[number_of_segments];
    sax_type region_upperv[number_of_segments];
    i=2;
    for (; i<number_of_segments; i++) 
    {
        region_lowerv[i] = (vv[i] <<  (c_m - c_cv[i]));
        region_upperv[i] = (~((int)MAXFLOAT << (c_m - c_cv[i])) | region_lowerv[i]);
    }

        float breakpoint_lowerv[number_of_segments] ;
        float breakpoint_upperv[number_of_segments] ;

    for (; i<number_of_segments; i++) 
    {
        if (region_lowerv[i] == 0) {
            breakpoint_lowerv[i] = min_val;
        }
        else
        {
            breakpoint_lowerv[i] = index->bins[i][  region_lowerv[i] - 1];
        }
        if (region_upperv[i] == max_cardinality - 1) {
            breakpoint_upperv[i] = max_val;
        }
        else
        {
            breakpoint_upperv[i] = index->bins[i][ region_upperv[i]];
        }

    }
ts_type valuev[number_of_segments];
    for (; i<number_of_segments; i++) 
    {
        if (breakpoint_lowerv[i] > fft[i]) {

            valuev[i] = breakpoint_lowerv[i] - fft[i];
        }
        else if(breakpoint_upperv[i] < fft[i]) {
            valuev[i] = (fft[i] - breakpoint_upperv[i]);
        }
    }

    for (; i<number_of_segments; i++) 
    {
        distance += 2*valuev[i]*valuev[i];
    }

    return distance;
}

ts_type minidist_fft_to_isax_raw_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf)
{
   
    int region_upper[16],region_lower[16];
    float distancef[16];
    int offset = 0;
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;

    __m256i vectorsignbit = _mm256_set1_epi32 (0xffffffff);



        //__m256i c_cv_0 = _mm256_set_epi32 ( sax_cardinalities[7] , sax_cardinalities[6] ,sax_cardinalities[5] ,sax_cardinalities[4] , sax_cardinalities[3] ,sax_cardinalities[2] ,sax_cardinalities[1],sax_cardinalities[0]);  
        //__m256i c_cv_1 = _mm256_set_epi32 ( sax_cardinalities[15], sax_cardinalities[14],sax_cardinalities[13],sax_cardinalities[12], sax_cardinalities[11],sax_cardinalities[10],sax_cardinalities[9],sax_cardinalities[8]);
        __m128i sax_cardinalitiesv8 = _mm_lddqu_si128 ((const void*)sax_cardinalities);
        __m256i sax_cardinalitiesv16= _mm256_cvtepu8_epi16 (sax_cardinalitiesv8);
        __m128i sax_cardinalitiesv16_0=_mm256_extractf128_si256 (sax_cardinalitiesv16,0);
        __m128i sax_cardinalitiesv16_1=_mm256_extractf128_si256 (sax_cardinalitiesv16,1);
        __m256i c_cv_0 = _mm256_cvtepu16_epi32 (sax_cardinalitiesv16_0);
        __m256i c_cv_1 = _mm256_cvtepu16_epi32 (sax_cardinalitiesv16_1);

        //__m256i v_0    = _mm256_set_epi32 (sax[7],sax[6],sax[5],sax[4],sax[3],sax[2],sax[1],sax[0]);
        //__m256i v_1    = _mm256_set_epi32 (sax[15],sax[14],sax[13],sax[12],sax[11],sax[10],sax[9],sax[8]);
        __m128i saxv8= _mm_lddqu_si128 ((const void*)sax);
        __m256i saxv16= _mm256_cvtepu8_epi16 (saxv8);
        __m128i saxv16_0 =_mm256_extractf128_si256 (saxv16,0);
        __m128i saxv16_1=_mm256_extractf128_si256 (saxv16,1);
        __m256i v_0= _mm256_cvtepu16_epi32 (saxv16_0);
        __m256i v_1 = _mm256_cvtepu16_epi32 (saxv16_1);


        __m256i c_m    = _mm256_set1_epi32 (max_bit_cardinality);
        __m256i cm_ccv_0 = _mm256_sub_epi32 (c_m, c_cv_0);
        __m256i cm_ccv_1 = _mm256_sub_epi32 (c_m, c_cv_1);

        //__m256i _mm256_set_epi32 (int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0)    
        //  __m256i _mm256_set1_epi32 (int a)
        __m256i region_lowerv_0 = _mm256_srlv_epi32 (v_0, cm_ccv_0);
        __m256i region_lowerv_1 = _mm256_srlv_epi32 (v_1, cm_ccv_1);
        region_lowerv_0 =  _mm256_sllv_epi32 (region_lowerv_0, cm_ccv_0);
        region_lowerv_1 =  _mm256_sllv_epi32 (region_lowerv_1, cm_ccv_1);

        
        __m256i v1= _mm256_andnot_si256 (_mm256_setzero_si256 (), vectorsignbit);

        __m256i region_upperv_0 = _mm256_sllv_epi32 (v1,cm_ccv_0);
        __m256i region_upperv_1 = _mm256_sllv_epi32 (v1,cm_ccv_1);
        region_upperv_0 = _mm256_andnot_si256 (region_upperv_0, vectorsignbit);
        region_upperv_1 = _mm256_andnot_si256 (region_upperv_1, vectorsignbit);

        region_upperv_0 = _mm256_or_si256 (region_upperv_0, region_lowerv_0);
        region_upperv_1 = _mm256_or_si256 (region_upperv_1, region_lowerv_1);

        _mm256_storeu_si256 ((void*)&(region_lower[0]),region_lowerv_0);
        _mm256_storeu_si256 ((void*)&(region_lower[8]),region_lowerv_1);
        _mm256_storeu_si256 ((void*)&(region_upper[0]),region_upperv_0);
        _mm256_storeu_si256 ((void*)&(region_upper[8]),region_upperv_1);

        
    //lower

        __m256i lower_juge_zerov_0 = _mm256_cmpeq_epi32 (region_lowerv_0, _mm256_setzero_si256 ());
        __m256i lower_juge_zerov_1 = _mm256_cmpeq_epi32 (region_lowerv_1, _mm256_setzero_si256 ());

        __m256i lower_juge_nzerov_0 = _mm256_andnot_si256 (lower_juge_zerov_0, vectorsignbit);
        __m256i lower_juge_nzerov_1 = _mm256_andnot_si256 (lower_juge_zerov_1, vectorsignbit);

        __m256 minvalv = _mm256_set1_ps (MINVAL);



        
        //__m256 lsax_breakpoints_shiftv_0 _mm256_i32gather_ps (sax_breakpoints, __m256i vindex, const int scale)
        __m256 lsax_breakpoints_shiftv_0= _mm256_set_ps (index->bins[7][ region_lower[7]-1],
                                                        index->bins[6][region_lower[6] - 1], 
                                                        index->bins[5][region_lower[5] - 1],
                                                        index->bins[4][region_lower[4] - 1],
                                                        index->bins[3][region_lower[3] - 1],
                                                        index->bins[2][region_lower[2] - 1],
                                                        index->bins[1][region_lower[1] - 1],
                                                        index->bins[0][region_lower[0] - 1]);
        __m256 lsax_breakpoints_shiftv_1= _mm256_set_ps (index->bins[15][region_lower[15] - 1],
                                                        index->bins[14][region_lower[14] - 1], 
                                                        index->bins[13][ region_lower[13] - 1],
                                                        index->bins[12][region_lower[12] - 1],
                                                        index->bins[11][region_lower[11] - 1],
                                                        index->bins[10][region_lower[10] - 1],
                                                        index->bins[9][region_lower[9] - 1],
                                                        index->bins[8][region_lower[8] - 1]);


        __m256 breakpoint_lowerv_0 = (__m256)_mm256_or_si256 (_mm256_and_si256(lower_juge_zerov_0,(__m256i)minvalv),_mm256_and_si256(lower_juge_nzerov_0,(__m256i)lsax_breakpoints_shiftv_0));
        __m256 breakpoint_lowerv_1 = (__m256)_mm256_or_si256 (_mm256_and_si256(lower_juge_zerov_1,(__m256i)minvalv),_mm256_and_si256(lower_juge_nzerov_1,(__m256i)lsax_breakpoints_shiftv_1));

    //uper
        __m256 usax_breakpoints_shiftv_0= _mm256_set_ps (index->bins[7][region_upper[7]],
                                                        index->bins[6][region_upper[6]], 
                                                        index->bins[5][region_upper[5]],
                                                        index->bins[4][region_upper[4]],
                                                        index->bins[3][region_upper[3]],
                                                        index->bins[2][region_upper[2]],
                                                        index->bins[1][region_upper[1]],
                                                        index->bins[0][region_upper[0]]);
        __m256 usax_breakpoints_shiftv_1= _mm256_set_ps (index->bins[15][region_upper[15]],
                                                        index->bins[14][region_upper[14]], 
                                                        index->bins[13][region_upper[13]],
                                                        index->bins[12][region_upper[12]],
                                                        index->bins[11][region_upper[11]],
                                                        index->bins[10][region_upper[10]],
                                                        index->bins[9][region_upper[9]],
                                                        index->bins[8][region_upper[8]]);


        __m256i upper_juge_maxv_0 = _mm256_cmpeq_epi32 (region_upperv_0,  _mm256_set1_epi32 (max_bit_cardinality - 1));
        __m256i upper_juge_maxv_1 = _mm256_cmpeq_epi32 (region_upperv_1,  _mm256_set1_epi32 (max_bit_cardinality - 1));

        __m256i upper_juge_nmaxv_0 = _mm256_andnot_si256 (upper_juge_maxv_0, vectorsignbit);
        __m256i upper_juge_nmaxv_1 = _mm256_andnot_si256 (upper_juge_maxv_1, vectorsignbit);

        __m256 breakpoint_upperv_0 = (__m256)_mm256_or_si256 (_mm256_and_si256(upper_juge_maxv_0,(__m256i)_mm256_set1_ps (MAXVAL)),_mm256_and_si256(upper_juge_nmaxv_0,(__m256i)usax_breakpoints_shiftv_0));
        __m256 breakpoint_upperv_1 = (__m256)_mm256_or_si256 (_mm256_and_si256(upper_juge_maxv_1,(__m256i)_mm256_set1_ps (MAXVAL)),_mm256_and_si256(upper_juge_nmaxv_1,(__m256i)usax_breakpoints_shiftv_1));







    //dis
            __m256 paav_0,paav_1;



            paav_0 =_mm256_loadu_ps (fft);
            paav_1 =_mm256_loadu_ps (&(fft[8]));
    


            __m256 dis_juge_upv_0=_mm256_cmp_ps (breakpoint_lowerv_0, paav_0, _CMP_GT_OS);
            __m256 dis_juge_upv_1=_mm256_cmp_ps (breakpoint_lowerv_1, paav_1, _CMP_GT_OS);




            __m256 dis_juge_lov_0=(__m256)_mm256_and_si256 ((__m256i)_mm256_cmp_ps (breakpoint_lowerv_0, paav_0, _CMP_NGT_US),(__m256i)_mm256_cmp_ps (breakpoint_upperv_0, paav_0, _CMP_LT_OS))  ;
            __m256 dis_juge_lov_1=(__m256)_mm256_and_si256 ((__m256i)_mm256_cmp_ps (breakpoint_lowerv_1, paav_1, _CMP_NGT_US),(__m256i)_mm256_cmp_ps (breakpoint_upperv_1, paav_1, _CMP_LT_OS));

            __m256 dis_juge_elv_0=(__m256)_mm256_andnot_si256 (_mm256_or_si256 ((__m256i)dis_juge_upv_0, (__m256i)dis_juge_lov_0),vectorsignbit);
            __m256 dis_juge_elv_1=(__m256)_mm256_andnot_si256 (_mm256_or_si256 ((__m256i)dis_juge_upv_1, (__m256i)dis_juge_lov_1),vectorsignbit);

            __m256 dis_lowv_0 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_lowerv_0, paav_0),_mm256_sub_ps (breakpoint_lowerv_0, paav_0));
            __m256 dis_lowv_1 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_lowerv_1, paav_1),_mm256_sub_ps (breakpoint_lowerv_1, paav_1));
            __m256 dis_uppv_0 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_upperv_0, paav_0),_mm256_sub_ps (breakpoint_upperv_0, paav_0));
            __m256 dis_uppv_1 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_upperv_1, paav_1),_mm256_sub_ps (breakpoint_upperv_1, paav_1));


            __m256 distancev_0=(__m256)_mm256_or_si256(_mm256_or_si256(_mm256_and_si256((__m256i)dis_juge_upv_0,(__m256i)dis_lowv_0),_mm256_and_si256((__m256i)dis_juge_lov_0,(__m256i)dis_uppv_0)),_mm256_and_si256((__m256i)dis_juge_elv_0,(__m256i)_mm256_set1_ps (0.0)));
            __m256 distancev_1=(__m256)_mm256_or_si256(_mm256_or_si256(_mm256_and_si256((__m256i)dis_juge_upv_1,(__m256i)dis_lowv_1),_mm256_and_si256((__m256i)dis_juge_lov_1,(__m256i)dis_uppv_1)),_mm256_and_si256((__m256i)dis_juge_elv_1,(__m256i)_mm256_set1_ps (0.0)));

            __m256 distancev = _mm256_add_ps (distancev_0, distancev_1);
            __m256 distancev2 = _mm256_hadd_ps (distancev, distancev);
            __m256 distancevf = _mm256_hadd_ps (distancev2, distancev2);


            _mm256_storeu_ps (distancef ,distancevf);
        //_mm256_storeu_ps (&checkvalue[8] ,distancev_1);


            return (distancef[0]+distancef[4])*2 ;
            
}

ts_type minidist_fft_to_isax_rawa_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf)
{
   
    int region_upper[16],region_lower[16];
    float distancef[8];
    int offset = 0;
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;

    __m256i vectorsignbit = _mm256_set1_epi32 (0xffffffff);



        //__m256i c_cv_0 = _mm256_set_epi32 ( sax_cardinalities[7] , sax_cardinalities[6] ,sax_cardinalities[5] ,sax_cardinalities[4] , sax_cardinalities[3] ,sax_cardinalities[2] ,sax_cardinalities[1],sax_cardinalities[0]);  
        //__m256i c_cv_1 = _mm256_set_epi32 ( sax_cardinalities[15], sax_cardinalities[14],sax_cardinalities[13],sax_cardinalities[12], sax_cardinalities[11],sax_cardinalities[10],sax_cardinalities[9],sax_cardinalities[8]);
        __m128i sax_cardinalitiesv8 = _mm_lddqu_si128 ((const void*)sax_cardinalities);
        __m256i sax_cardinalitiesv16= _mm256_cvtepu8_epi16 (sax_cardinalitiesv8);
        __m128i sax_cardinalitiesv16_0=_mm256_extractf128_si256 (sax_cardinalitiesv16,0);
        __m128i sax_cardinalitiesv16_1=_mm256_extractf128_si256 (sax_cardinalitiesv16,1);
        __m256i c_cv_0 = _mm256_cvtepu16_epi32 (sax_cardinalitiesv16_0);
        __m256i c_cv_1 = _mm256_cvtepu16_epi32 (sax_cardinalitiesv16_1);

        //__m256i v_0    = _mm256_set_epi32 (sax[7],sax[6],sax[5],sax[4],sax[3],sax[2],sax[1],sax[0]);
        //__m256i v_1    = _mm256_set_epi32 (sax[15],sax[14],sax[13],sax[12],sax[11],sax[10],sax[9],sax[8]);
        __m128i saxv8= _mm_lddqu_si128 ((const void*)sax);
        __m256i saxv16= _mm256_cvtepu8_epi16 (saxv8);
        __m128i saxv16_0 =_mm256_extractf128_si256 (saxv16,0);
        __m128i saxv16_1=_mm256_extractf128_si256 (saxv16,1);
        __m256i v_0= _mm256_cvtepu16_epi32 (saxv16_0);
        __m256i v_1 = _mm256_cvtepu16_epi32 (saxv16_1);


        __m256i c_m    = _mm256_set1_epi32 (max_bit_cardinality);
        __m256i cm_ccv_0 = _mm256_sub_epi32 (c_m, c_cv_0);
        __m256i cm_ccv_1 = _mm256_sub_epi32 (c_m, c_cv_1);

        //__m256i _mm256_set_epi32 (int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0)    
        //  __m256i _mm256_set1_epi32 (int a)
        __m256i region_lowerv_0 = _mm256_srlv_epi32 (v_0, cm_ccv_0);
        __m256i region_lowerv_1 = _mm256_srlv_epi32 (v_1, cm_ccv_1);
        region_lowerv_0 =  _mm256_sllv_epi32 (region_lowerv_0, cm_ccv_0);
        region_lowerv_1 =  _mm256_sllv_epi32 (region_lowerv_1, cm_ccv_1);

        
        __m256i v1= _mm256_andnot_si256 (_mm256_setzero_si256 (), vectorsignbit);

        __m256i region_upperv_0 = _mm256_sllv_epi32 (v1,cm_ccv_0);
        __m256i region_upperv_1 = _mm256_sllv_epi32 (v1,cm_ccv_1);
        region_upperv_0 = _mm256_andnot_si256 (region_upperv_0, vectorsignbit);
        region_upperv_1 = _mm256_andnot_si256 (region_upperv_1, vectorsignbit);

        region_upperv_0 = _mm256_or_si256 (region_upperv_0, region_lowerv_0);
        region_upperv_1 = _mm256_or_si256 (region_upperv_1, region_lowerv_1);

        //_mm256_storeu_si256 ((void*)&(region_lower[0]),region_lowerv_0);
        //_mm256_storeu_si256 ((void*)&(region_lower[8]),region_lowerv_1);
        //_mm256_storeu_si256 ((void*)&(region_upper[0]),region_upperv_0);
        //_mm256_storeu_si256 ((void*)&(region_upper[8]),region_upperv_1);

        
    //lower

        __m256i lower_juge_zerov_0 = _mm256_cmpeq_epi32 (region_lowerv_0, _mm256_setzero_si256 ());
        __m256i lower_juge_zerov_1 = _mm256_cmpeq_epi32 (region_lowerv_1, _mm256_setzero_si256 ());

        __m256i lower_juge_nzerov_0 = _mm256_andnot_si256 (lower_juge_zerov_0, vectorsignbit);
        __m256i lower_juge_nzerov_1 = _mm256_andnot_si256 (lower_juge_zerov_1, vectorsignbit);

        __m256 minvalv = _mm256_set1_ps (MINVAL);
        __m256i bitsizev = _mm256_set1_epi16 ((short)index->settings->sax_alphabet_cardinality-1);
        __m256i bit1v = _mm256_set1_epi32(1);
        //__m256i bit8v = _mm256_set1_epi32(8);
        __m256i offsetvs=_mm256_set_epi16(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
        //__m256i offsetv1=_mm256_set_epi32();
        //__m256i offsetv1=_mm256_add_epi32(offsetv0,bit8v);

       __m256i vsssssoffsetvs2= _mm256_mullo_epi16 (bitsizev, offsetvs);
        __m128i offsetv0s= _mm256_extractf128_si256 (vsssssoffsetvs2, 0);
        __m128i offsetv1s= _mm256_extractf128_si256 (vsssssoffsetvs2, 1);
        __m256i  offsetv0=   _mm256_cvtepu16_epi32(offsetv0s);
        __m256i  offsetv1=   _mm256_cvtepu16_epi32(offsetv1s);
        //__m256i offsetvv0=_mm256_mul_epi32 (offsetv0, bitsizev);
        //__m256i offsetvv1=_mm256_mul_epi32 (offsetv1, bitsizev);

        
        __m256i region_lowerbinv0=_mm256_add_epi32 (offsetv0, region_lowerv_0);
        __m256i region_lowerbinv1=_mm256_add_epi32 (offsetv1, region_lowerv_1);
        region_lowerbinv0=_mm256_sub_epi32 (region_lowerbinv0,bit1v);
        region_lowerbinv1=_mm256_sub_epi32 (region_lowerbinv1,bit1v);
        //__m256 lsax_breakpoints_shiftv_0 _mm256_i32gather_ps (sax_breakpoints, __m256i vindex, const int scale)
        __m256 lsax_breakpoints_shiftv_0 =_mm256_i32gather_ps (index->binsv,region_lowerbinv0, 4);
        /*__m256 lsax_breakpoints_shiftv_0= _mm256_set_ps (index->binsv[region_lower[7]+7*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[6]+6*(index->settings->sax_alphabet_cardinality-1) -1], 
                                                        index->binsv[region_lower[5]+5*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[4]+4*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[3]+3*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[2]+2*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[1]+1*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[0] - 1]);
*/
        __m256 lsax_breakpoints_shiftv_1 =_mm256_i32gather_ps (index->binsv,region_lowerbinv1, 4);
        /*__m256 lsax_breakpoints_shiftv_1= _mm256_set_ps (index->binsv[region_lower[15]+15*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[14] +14*(index->settings->sax_alphabet_cardinality-1)-1], 
                                                        index->binsv[region_lower[13]+13*(index->settings->sax_alphabet_cardinality-1)-1],
                                                        index->binsv[region_lower[12] +12*(index->settings->sax_alphabet_cardinality-1)-1],
                                                        index->binsv[region_lower[11]+11*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[10]+10*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[9] +9*(index->settings->sax_alphabet_cardinality-1)-1],
                                                        index->binsv[region_lower[8] +8*(index->settings->sax_alphabet_cardinality-1)-1]);

*/
        __m256 breakpoint_lowerv_0 = (__m256)_mm256_or_si256 (_mm256_and_si256(lower_juge_zerov_0,(__m256i)minvalv),_mm256_and_si256(lower_juge_nzerov_0,(__m256i)lsax_breakpoints_shiftv_0));
        __m256 breakpoint_lowerv_1 = (__m256)_mm256_or_si256 (_mm256_and_si256(lower_juge_zerov_1,(__m256i)minvalv),_mm256_and_si256(lower_juge_nzerov_1,(__m256i)lsax_breakpoints_shiftv_1));

    //uper        
        __m256i region_upperbinv0=_mm256_add_epi32 (offsetv0, region_upperv_0);
        __m256i region_upperbinv1=_mm256_add_epi32 (offsetv1, region_upperv_1);
        __m256 usax_breakpoints_shiftv_0 =_mm256_i32gather_ps (index->binsv,region_upperbinv0, 4);
        /*__m256 usax_breakpoints_shiftv_0= _mm256_set_ps (index->binsv[region_upper[7]+7*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[6]+6*(index->settings->sax_alphabet_cardinality-1)], 
                                                        index->binsv[region_upper[5]+5*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[4]+4*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[3]+3*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[2]+2*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[1]+1*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[0]+0*(index->settings->sax_alphabet_cardinality-1)]);
        */__m256 usax_breakpoints_shiftv_1 =_mm256_i32gather_ps (index->binsv,region_upperbinv1, 4);
        /*__m256 usax_breakpoints_shiftv_1= _mm256_set_ps (index->binsv[region_upper[15]+15*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[14]+14*(index->settings->sax_alphabet_cardinality-1)], 
                                                        index->binsv[region_upper[13]+13*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[12]+12*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[11]+11*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[10]+10*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[9]+9*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[8]+8*(index->settings->sax_alphabet_cardinality-1)]);
*/

        __m256i upper_juge_maxv_0 = _mm256_cmpeq_epi32 (region_upperv_0,  _mm256_set1_epi32 (max_bit_cardinality - 1));
        __m256i upper_juge_maxv_1 = _mm256_cmpeq_epi32 (region_upperv_1,  _mm256_set1_epi32 (max_bit_cardinality - 1));

        __m256i upper_juge_nmaxv_0 = _mm256_andnot_si256 (upper_juge_maxv_0, vectorsignbit);
        __m256i upper_juge_nmaxv_1 = _mm256_andnot_si256 (upper_juge_maxv_1, vectorsignbit);

        __m256 breakpoint_upperv_0 = (__m256)_mm256_or_si256 (_mm256_and_si256(upper_juge_maxv_0,(__m256i)_mm256_set1_ps (MAXVAL)),_mm256_and_si256(upper_juge_nmaxv_0,(__m256i)usax_breakpoints_shiftv_0));
        __m256 breakpoint_upperv_1 = (__m256)_mm256_or_si256 (_mm256_and_si256(upper_juge_maxv_1,(__m256i)_mm256_set1_ps (MAXVAL)),_mm256_and_si256(upper_juge_nmaxv_1,(__m256i)usax_breakpoints_shiftv_1));







    //dis
            __m256 paav_0,paav_1;



            paav_0 =_mm256_loadu_ps (fft);
            paav_1 =_mm256_loadu_ps (&(fft[8]));
    


            __m256 dis_juge_upv_0=_mm256_cmp_ps (breakpoint_lowerv_0, paav_0, _CMP_GT_OS);
            __m256 dis_juge_upv_1=_mm256_cmp_ps (breakpoint_lowerv_1, paav_1, _CMP_GT_OS);




            __m256 dis_juge_lov_0=(__m256)_mm256_and_si256 ((__m256i)_mm256_cmp_ps (breakpoint_lowerv_0, paav_0, _CMP_NGT_US),(__m256i)_mm256_cmp_ps (breakpoint_upperv_0, paav_0, _CMP_LT_OS))  ;
            __m256 dis_juge_lov_1=(__m256)_mm256_and_si256 ((__m256i)_mm256_cmp_ps (breakpoint_lowerv_1, paav_1, _CMP_NGT_US),(__m256i)_mm256_cmp_ps (breakpoint_upperv_1, paav_1, _CMP_LT_OS));

            __m256 dis_juge_elv_0=(__m256)_mm256_andnot_si256 (_mm256_or_si256 ((__m256i)dis_juge_upv_0, (__m256i)dis_juge_lov_0),vectorsignbit);
            __m256 dis_juge_elv_1=(__m256)_mm256_andnot_si256 (_mm256_or_si256 ((__m256i)dis_juge_upv_1, (__m256i)dis_juge_lov_1),vectorsignbit);

            __m256 dis_lowv_0 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_lowerv_0, paav_0),_mm256_sub_ps (breakpoint_lowerv_0, paav_0));
            __m256 dis_lowv_1 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_lowerv_1, paav_1),_mm256_sub_ps (breakpoint_lowerv_1, paav_1));
            __m256 dis_uppv_0 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_upperv_0, paav_0),_mm256_sub_ps (breakpoint_upperv_0, paav_0));
            __m256 dis_uppv_1 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_upperv_1, paav_1),_mm256_sub_ps (breakpoint_upperv_1, paav_1));


            __m256 distancev_0=(__m256)_mm256_or_si256(_mm256_or_si256(_mm256_and_si256((__m256i)dis_juge_upv_0,(__m256i)dis_lowv_0),_mm256_and_si256((__m256i)dis_juge_lov_0,(__m256i)dis_uppv_0)),_mm256_and_si256((__m256i)dis_juge_elv_0,(__m256i)_mm256_set1_ps (0.0)));
            __m256 distancev_1=(__m256)_mm256_or_si256(_mm256_or_si256(_mm256_and_si256((__m256i)dis_juge_upv_1,(__m256i)dis_lowv_1),_mm256_and_si256((__m256i)dis_juge_lov_1,(__m256i)dis_uppv_1)),_mm256_and_si256((__m256i)dis_juge_elv_1,(__m256i)_mm256_set1_ps (0.0)));

            __m256 distancev = _mm256_add_ps (distancev_0, distancev_1);
            __m256 distancev2 = _mm256_hadd_ps (distancev, distancev);
            __m256 distancevf = _mm256_hadd_ps (distancev2, distancev2);


            _mm256_storeu_ps (distancef ,distancevf);
        //_mm256_storeu_ps (&checkvalue[8] ,distancev_1);
                    //_mm256_storeu_si256 ((void*)&(offsetvvvv[0]),offsetv0);
       // _mm256_storeu_si256 ((void*)&(offsetvvvv[8]),offsetv1);
        //for (int i = 0; i < 16; i++)
        //{
          //  printf("the number is victor [%d] %d\n",i,offsetvvvv[i]);
           // printf("the number is [%d] %d\n",i,i*(index->settings->sax_alphabet_cardinality-1));

       // }
       //printf("distance is %f!!!!!!!\n",(distancef[0]+distancef[4])*2);
        //sleep(5);

            return (distancef[0]+distancef[4])*2 ;
            
}


ts_type minidist_fft_to_isax_rawe_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf)
{
   
    int region_upper[16],region_lower[16];
    float distancef[8],distancef2[8];
    int offset = 0;
    sax_type max_bit_cardinality = index->settings->sax_bit_cardinality;

    __m256i vectorsignbit = _mm256_set1_epi32 (0xffffffff);



        //__m256i c_cv_0 = _mm256_set_epi32 ( sax_cardinalities[7] , sax_cardinalities[6] ,sax_cardinalities[5] ,sax_cardinalities[4] , sax_cardinalities[3] ,sax_cardinalities[2] ,sax_cardinalities[1],sax_cardinalities[0]);  
        //__m256i c_cv_1 = _mm256_set_epi32 ( sax_cardinalities[15], sax_cardinalities[14],sax_cardinalities[13],sax_cardinalities[12], sax_cardinalities[11],sax_cardinalities[10],sax_cardinalities[9],sax_cardinalities[8]);
        __m128i sax_cardinalitiesv8 = _mm_lddqu_si128 ((const void*)sax_cardinalities);
        __m256i sax_cardinalitiesv16= _mm256_cvtepu8_epi16 (sax_cardinalitiesv8);
        __m128i sax_cardinalitiesv16_0=_mm256_extractf128_si256 (sax_cardinalitiesv16,0);
        __m256i c_cv_0 = _mm256_cvtepu16_epi32 (sax_cardinalitiesv16_0);
 
        //__m256i v_0    = _mm256_set_epi32 (sax[7],sax[6],sax[5],sax[4],sax[3],sax[2],sax[1],sax[0]);
        //__m256i v_1    = _mm256_set_epi32 (sax[15],sax[14],sax[13],sax[12],sax[11],sax[10],sax[9],sax[8]);
        __m128i saxv8= _mm_lddqu_si128 ((const void*)sax);
        __m256i saxv16= _mm256_cvtepu8_epi16 (saxv8);
        __m128i saxv16_0 =_mm256_extractf128_si256 (saxv16,0);

        __m256i v_0= _mm256_cvtepu16_epi32 (saxv16_0);



        __m256i c_m    = _mm256_set1_epi32 (max_bit_cardinality);
        __m256i cm_ccv_0 = _mm256_sub_epi32 (c_m, c_cv_0);


        //__m256i _mm256_set_epi32 (int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0)    
        //  __m256i _mm256_set1_epi32 (int a)
        __m256i region_lowerv_0 = _mm256_srlv_epi32 (v_0, cm_ccv_0);

        region_lowerv_0 =  _mm256_sllv_epi32 (region_lowerv_0, cm_ccv_0);


        
        __m256i v1= _mm256_andnot_si256 (_mm256_setzero_si256 (), vectorsignbit);

        __m256i region_upperv_0 = _mm256_sllv_epi32 (v1,cm_ccv_0);

        region_upperv_0 = _mm256_andnot_si256 (region_upperv_0, vectorsignbit);


        region_upperv_0 = _mm256_or_si256 (region_upperv_0, region_lowerv_0);


        //_mm256_storeu_si256 ((void*)&(region_lower[0]),region_lowerv_0);
        //_mm256_storeu_si256 ((void*)&(region_lower[8]),region_lowerv_1);
        //_mm256_storeu_si256 ((void*)&(region_upper[0]),region_upperv_0);
        //_mm256_storeu_si256 ((void*)&(region_upper[8]),region_upperv_1);

        
    //lower

        __m256i lower_juge_zerov_0 = _mm256_cmpeq_epi32 (region_lowerv_0, _mm256_setzero_si256 ());


        __m256i lower_juge_nzerov_0 = _mm256_andnot_si256 (lower_juge_zerov_0, vectorsignbit);


        __m256 minvalv = _mm256_set1_ps (MINVAL);
        __m256i bitsizev = _mm256_set1_epi16 ((short)index->settings->sax_alphabet_cardinality-1);
        __m256i bit1v = _mm256_set1_epi32(1);
        //__m256i bit8v = _mm256_set1_epi32(8);
        __m256i offsetvs=_mm256_set_epi16(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0);
        //__m256i offsetv1=_mm256_set_epi32();
        //__m256i offsetv1=_mm256_add_epi32(offsetv0,bit8v);

       __m256i vsssssoffsetvs2= _mm256_mullo_epi16 (bitsizev, offsetvs);
        __m128i offsetv0s= _mm256_extractf128_si256 (vsssssoffsetvs2, 0);

        __m128i offsetv1s= _mm256_extractf128_si256 (vsssssoffsetvs2, 1);
        __m256i  offsetv0=   _mm256_cvtepu16_epi32(offsetv0s);

        //__m256i offsetvv0=_mm256_mul_epi32 (offsetv0, bitsizev);
        //__m256i offsetvv1=_mm256_mul_epi32 (offsetv1, bitsizev);

        
        __m256i region_lowerbinv0=_mm256_add_epi32 (offsetv0, region_lowerv_0);

        region_lowerbinv0=_mm256_sub_epi32 (region_lowerbinv0,bit1v);

        //__m256 lsax_breakpoints_shiftv_0 _mm256_i32gather_ps (sax_breakpoints, __m256i vindex, const int scale)
        __m256 lsax_breakpoints_shiftv_0 =_mm256_i32gather_ps (index->binsv,region_lowerbinv0, 4);
        /*__m256 lsax_breakpoints_shiftv_0= _mm256_set_ps (index->binsv[region_lower[7]+7*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[6]+6*(index->settings->sax_alphabet_cardinality-1) -1], 
                                                        index->binsv[region_lower[5]+5*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[4]+4*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[3]+3*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[2]+2*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[1]+1*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[0] - 1]);
*/

        /*__m256 lsax_breakpoints_shiftv_1= _mm256_set_ps (index->binsv[region_lower[15]+15*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[14] +14*(index->settings->sax_alphabet_cardinality-1)-1], 
                                                        index->binsv[region_lower[13]+13*(index->settings->sax_alphabet_cardinality-1)-1],
                                                        index->binsv[region_lower[12] +12*(index->settings->sax_alphabet_cardinality-1)-1],
                                                        index->binsv[region_lower[11]+11*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[10]+10*(index->settings->sax_alphabet_cardinality-1) -1],
                                                        index->binsv[region_lower[9] +9*(index->settings->sax_alphabet_cardinality-1)-1],
                                                        index->binsv[region_lower[8] +8*(index->settings->sax_alphabet_cardinality-1)-1]);

*/
        __m256 breakpoint_lowerv_0 = (__m256)_mm256_or_si256 (_mm256_and_si256(lower_juge_zerov_0,(__m256i)minvalv),_mm256_and_si256(lower_juge_nzerov_0,(__m256i)lsax_breakpoints_shiftv_0));

    //uper        
        __m256i region_upperbinv0=_mm256_add_epi32 (offsetv0, region_upperv_0);

        __m256 usax_breakpoints_shiftv_0 =_mm256_i32gather_ps (index->binsv,region_upperbinv0, 4);
        /*__m256 usax_breakpoints_shiftv_0= _mm256_set_ps (index->binsv[region_upper[7]+7*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[6]+6*(index->settings->sax_alphabet_cardinality-1)], 
                                                        index->binsv[region_upper[5]+5*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[4]+4*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[3]+3*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[2]+2*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[1]+1*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[0]+0*(index->settings->sax_alphabet_cardinality-1)]);
        */

        /*__m256 usax_breakpoints_shiftv_1= _mm256_set_ps (index->binsv[region_upper[15]+15*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[14]+14*(index->settings->sax_alphabet_cardinality-1)], 
                                                        index->binsv[region_upper[13]+13*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[12]+12*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[11]+11*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[10]+10*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[9]+9*(index->settings->sax_alphabet_cardinality-1)],
                                                        index->binsv[region_upper[8]+8*(index->settings->sax_alphabet_cardinality-1)]);
*/

        __m256i upper_juge_maxv_0 = _mm256_cmpeq_epi32 (region_upperv_0,  _mm256_set1_epi32 (max_bit_cardinality - 1));

        __m256i upper_juge_nmaxv_0 = _mm256_andnot_si256 (upper_juge_maxv_0, vectorsignbit);

        __m256 breakpoint_upperv_0 = (__m256)_mm256_or_si256 (_mm256_and_si256(upper_juge_maxv_0,(__m256i)_mm256_set1_ps (MAXVAL)),_mm256_and_si256(upper_juge_nmaxv_0,(__m256i)usax_breakpoints_shiftv_0));






    //dis
            __m256 paav_0,paav_1;



            paav_0 =_mm256_loadu_ps (fft);



            __m256 dis_juge_upv_0=_mm256_cmp_ps (breakpoint_lowerv_0, paav_0, _CMP_GT_OS);
                                        



            __m256 dis_juge_lov_0=(__m256)_mm256_and_si256 ((__m256i)_mm256_cmp_ps (breakpoint_lowerv_0, paav_0, _CMP_NGT_US),(__m256i)_mm256_cmp_ps (breakpoint_upperv_0, paav_0, _CMP_LT_OS))  ;

            __m256 dis_juge_elv_0=(__m256)_mm256_andnot_si256 (_mm256_or_si256 ((__m256i)dis_juge_upv_0, (__m256i)dis_juge_lov_0),vectorsignbit);

            __m256 dis_lowv_0 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_lowerv_0, paav_0),_mm256_sub_ps (breakpoint_lowerv_0, paav_0));
 __m256 dis_uppv_0 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_upperv_0, paav_0),_mm256_sub_ps (breakpoint_upperv_0, paav_0));
            


            __m256 distancev_0=(__m256)_mm256_or_si256(_mm256_or_si256(_mm256_and_si256((__m256i)dis_juge_upv_0,(__m256i)dis_lowv_0),_mm256_and_si256((__m256i)dis_juge_lov_0,(__m256i)dis_uppv_0)),_mm256_and_si256((__m256i)dis_juge_elv_0,(__m256i)_mm256_set1_ps (0.0)));

            __m256 distancev2 = _mm256_hadd_ps (distancev_0, distancev_0);
            __m256 distancevf = _mm256_hadd_ps (distancev2, distancev2);


            _mm256_storeu_ps (distancef ,distancevf);
        //_mm256_storeu_ps (&checkvalue[8] ,distancev_1);
                    //_mm256_storeu_si256 ((void*)&(offsetvvvv[0]),offsetv0);
       // _mm256_storeu_si256 ((void*)&(offsetvvvv[8]),offsetv1);
        //for (int i = 0; i < 16; i++)
        //{
          //  printf("the number is victor [%d] %d\n",i,offsetvvvv[i]);
           // printf("the number is [%d] %d\n",i,i*(index->settings->sax_alphabet_cardinality-1));

       // }
       //printf("distance is %f!!!!!!!\n",(distancef[0]+distancef[4])*2);
        //sleep(5);
        if ( (distancef[0]+distancef[4])*2 >bsf)
        {
            return (distancef[0]+distancef[4])*2 ;
        }
        

        __m128i sax_cardinalitiesv16_1=_mm256_extractf128_si256 (sax_cardinalitiesv16,1);
        __m256i c_cv_1 = _mm256_cvtepu16_epi32 (sax_cardinalitiesv16_1);
        __m128i saxv16_1=_mm256_extractf128_si256 (saxv16,1);
        __m256i v_1 = _mm256_cvtepu16_epi32 (saxv16_1);
        __m256i cm_ccv_1 = _mm256_sub_epi32 (c_m, c_cv_1);
        __m256i region_lowerv_1 = _mm256_srlv_epi32 (v_1, cm_ccv_1);
        region_lowerv_1 =  _mm256_sllv_epi32 (region_lowerv_1, cm_ccv_1);
        __m256i region_upperv_1 = _mm256_sllv_epi32 (v1,cm_ccv_1);
        region_upperv_1 = _mm256_andnot_si256 (region_upperv_1, vectorsignbit);
        region_upperv_1 = _mm256_or_si256 (region_upperv_1, region_lowerv_1);
        __m256i lower_juge_zerov_1 = _mm256_cmpeq_epi32 (region_lowerv_1, _mm256_setzero_si256 ());
        __m256i lower_juge_nzerov_1 = _mm256_andnot_si256 (lower_juge_zerov_1, vectorsignbit);
        __m256i  offsetv1=   _mm256_cvtepu16_epi32(offsetv1s);
        __m256i region_lowerbinv1=_mm256_add_epi32 (offsetv1, region_lowerv_1);
        region_lowerbinv1=_mm256_sub_epi32 (region_lowerbinv1,bit1v);
        __m256 lsax_breakpoints_shiftv_1 =_mm256_i32gather_ps (index->binsv,region_lowerbinv1, 4);
        __m256 breakpoint_lowerv_1 = (__m256)_mm256_or_si256 (_mm256_and_si256(lower_juge_zerov_1,(__m256i)minvalv),_mm256_and_si256(lower_juge_nzerov_1,(__m256i)lsax_breakpoints_shiftv_1));

        __m256i region_upperbinv1=_mm256_add_epi32 (offsetv1, region_upperv_1);
        __m256 usax_breakpoints_shiftv_1 =_mm256_i32gather_ps (index->binsv,region_upperbinv1, 4);
        __m256i upper_juge_maxv_1 = _mm256_cmpeq_epi32 (region_upperv_1,  _mm256_set1_epi32 (max_bit_cardinality - 1));
        __m256i upper_juge_nmaxv_1 = _mm256_andnot_si256 (upper_juge_maxv_1, vectorsignbit);
        __m256 breakpoint_upperv_1 = (__m256)_mm256_or_si256 (_mm256_and_si256(upper_juge_maxv_1,(__m256i)_mm256_set1_ps (MAXVAL)),_mm256_and_si256(upper_juge_nmaxv_1,(__m256i)usax_breakpoints_shiftv_1));
        paav_1 =_mm256_loadu_ps (&(fft[8]));    
        __m256 dis_juge_upv_1=_mm256_cmp_ps (breakpoint_lowerv_1, paav_1, _CMP_GT_OS);


            __m256 dis_juge_lov_1=(__m256)_mm256_and_si256 ((__m256i)_mm256_cmp_ps (breakpoint_lowerv_1, paav_1, _CMP_NGT_US),(__m256i)_mm256_cmp_ps (breakpoint_upperv_1, paav_1, _CMP_LT_OS));

            __m256 dis_juge_elv_1=(__m256)_mm256_andnot_si256 (_mm256_or_si256 ((__m256i)dis_juge_upv_1, (__m256i)dis_juge_lov_1),vectorsignbit);

            __m256 dis_lowv_1 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_lowerv_1, paav_1),_mm256_sub_ps (breakpoint_lowerv_1, paav_1));
           
            __m256 dis_uppv_1 =_mm256_mul_ps (_mm256_sub_ps (breakpoint_upperv_1, paav_1),_mm256_sub_ps (breakpoint_upperv_1, paav_1));

            __m256 distancev_1=(__m256)_mm256_or_si256(_mm256_or_si256(_mm256_and_si256((__m256i)dis_juge_upv_1,(__m256i)dis_lowv_1),_mm256_and_si256((__m256i)dis_juge_lov_1,(__m256i)dis_uppv_1)),_mm256_and_si256((__m256i)dis_juge_elv_1,(__m256i)_mm256_set1_ps (0.0)));



             distancev2 = _mm256_hadd_ps (distancev_1, distancev_1);
            distancevf = _mm256_hadd_ps (distancev2, distancev2);
            _mm256_storeu_ps (distancef2 ,distancevf);

            return (distancef[0]+distancef[4]+distancef2[0]+distancef2[4])*2 ;
            
}