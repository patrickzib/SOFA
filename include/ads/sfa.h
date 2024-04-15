#ifndef sfalib_h
#define sfalib_h

#include "../../config.h"
#include "../../globals.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include <fftw3.h>
#include <sys/types.h>

enum response sfa_bins_init(isax_index *index);
void sfa_free_bins(isax_index *index);

void sfa_set_bins(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread, int filetype_int, int apply_znorm);
void sfa_set_bins_coeff(isax_index *index, const char *ifilename, long int ts_num, int maxquerythread, int filetype_int, int apply_znorm);

ts_type** calculate_variance_coeff(isax_index *index, ts_type ** dft_mem_array);

void* set_bins_worker_dft(void *transferdata);
void* set_bins_worker_dft_coeff(void *transferdata);

void* order_divide_worker(void *transferdata);

void sfa_print_bins(isax_index *index);

void free_dft_memory(isax_index *index, ts_type **dft_mem_array);

int compare_ts_type (const void * a, const void * b);
int compare_var (const void *a, const void *b);
int compare_int (const void *a, const void *b);

ts_type minidist_fft_to_isax(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_fft_to_isax_raw(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_fft_to_isax_raw_autoSIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_fft_to_isax_raw_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_fft_to_isax_rawa_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf);
ts_type minidist_fft_to_isax_rawe_SIMD(isax_index *index, float *fft, sax_type *sax, sax_type *sax_cardinalities, float bsf);

long random_at_most(long max);

void sfa_printbin(unsigned long long n, int size);

void fft_print(ts_type *fft, int segments);

void sfa_print(sax_type *sax,int segments,  int cardinality);	

typedef struct bins_data_inmemory
{
	isax_index *index;
	long int start_number,stop_number;	
	ts_type ** dft_mem_array;
    const char *filename;
	int workernumber;
	long int records;
	long int records_offset;
	ts_type * ts;
	fftwf_complex *ts_out;
	fftwf_plan plan_forward;
	ts_type * transform;
    int filetype_int;
    int apply_znorm;
}bins_data_inmemory;

typedef struct variance_coeff_index
{
	double variance;
	int coeff_index;
}variance_coeff_index;


#endif
