//
//  dft.c
//  sfa C version
//
//  Created by Karima Echihabi on 18/12/2016
//  Copyright 2016 Paris Descartes University. All rights reserved.
//  
//

#include <stdio.h>
#include <stdlib.h>
#include "../../config.h"
#include "../../globals.h"
#include "ads/isax_index.h"
#include "math.h"

#include <fftw3.h>

#include "ads/dft.h"

//void fft_from_ts_chunk(struct sfa_trie * trie)


double real_ephi (double u, double m)
{
  return cos(2*M_PI*u/m);
}

double complex_ephi (double u, double m)
{
  return -sin(2*M_PI*u/m);
}

double complex_mul_real(double r1, double im1, double r2, double im2)
{
    return r1*r2 - im1*im2;
}

double complex_mul_imag(double r1, double im1, double r2, double im2)
{
    return r1*im2 + r2*im1;
}

void fft_from_ts(isax_index *index, ts_type *ts, fftwf_complex *ts_out, ts_type *transform, fftwf_plan plan_forward)
{
	//fprintf(stderr,"fft from ts\n");

    unsigned long ts_length = index->settings->timeseries_size;
    int transforms_size = index->settings->fft_size;
      
    fftwf_execute (plan_forward);
	   
    ts_out[0][1] = 0;

    //check if in the case of whole matching, we need to multiply by norm and -1 
    int j = 0;
	   
    //if normalized, ignore first coeff and start with offset 1
    int start_offset = index->settings->is_norm?1:0;

    for(int k = start_offset; k< transforms_size/2+start_offset; ++k)
    {
	    transform[j] = ts_out[k][0];
	    transform[j+1] =  ts_out[k][1];
	    j +=2;
    }

    //normalizing fft result in frequency domain
    //TODO check is normalization with sign is right???
    int sign = 1;
    ts_type norm_factor = index->norm_factor;
       
    for(int i = 0; i< transforms_size; ++i)
    {
  	    transform[i] *= norm_factor*sign;
	    sign *= -1;
    }
    
    return;
}

/*
  The current transform is pointed to by dft_mem_array
*/

void sfa_from_fft(isax_index * index, ts_type * cur_transform, unsigned char * cur_sfa_word)
{
	//fprintf(stderr,"sfa from fft\n");

    unsigned long ts_length = index->settings->timeseries_size;
    int word_length = index->settings->fft_size;
      
    for (int k = 0; k < word_length; ++k)
    {
        unsigned int c;
        for (c = 0; c < index->settings->sax_alphabet_cardinality-1; c++)
        {
            if (cur_transform[k] < index->bins[k][c])
	        {
                break;
            }
        }
        cur_sfa_word[k] = (unsigned char) (c);
    }
}

enum response sfa_from_ts(isax_index *index, ts_type *ts_in, sax_type *sax_out, fftwf_complex *ts_out, ts_type *transform, fftwf_plan plan_forward)
{
    /*
    fprintf(stderr,"IDX TS word: "); 
    for(int i=0; i<index->settings->timeseries_size; ++i)
    {
        fprintf(stderr,"-%.3f", ts_in[i]); 
    }*/

    fft_from_ts(index, ts_in, ts_out, transform, plan_forward);

    ts_type * cur_coeff_line = calloc(index->settings->fft_size, sizeof(ts_type));
    
    for(int i=0; i<index->settings->fft_size; ++i)
    {
        cur_coeff_line[i] = (ts_type) roundf(transform[i]*100.0)/100.0;
    }
  
    /*
    fprintf(stderr,"IDX FFT word: "); 
    for(int i=0; i<index->settings->fft_size; ++i)
    {
        fprintf(stderr,"-%.3f", transform[i]); 
    }
    fprintf(stderr,"-\n");*/

    sfa_from_fft(index, cur_coeff_line, sax_out);
    /*
    fprintf(stderr,"SFA word: "); 
    for(int i=0; i<index->settings->fft_size; ++i)
    {
        fprintf(stderr,"-%u", sax_out[i]); 
    }
    fprintf(stderr,"-\n"); */

    if(sax_out != NULL) return SUCCESS;
    else
    {
        fprintf(stderr, "SFA error");
    }
}


