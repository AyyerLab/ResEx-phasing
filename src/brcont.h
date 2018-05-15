#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_rng.h>
#include <omp.h>
#include <sys/stat.h>
#include "fft.h"
#include "input.h"
#include "volume.h"

// Algorithm parameters
struct algorithm_data {
	long size, vol, num_vox ;
	char output_prefix[1024] ;
	
	// Other data structs
	struct volume_data *volume ;
	struct input_data *input ;
	struct fft_data *fft ;
	struct rotation *quat ;
	
	// Temporary arrays
	float *exp_mag, *iterate ;
	float *p1, *p2, *r1, *r2;
	float *average_p1, *average_p2 ;
	
	// Algorithms per iterate
	int num_iter, num_avg_iter ;
	char (*algorithms)[8], (*avg_algorithms)[8] ;
	
	// Parameters
	float beta ;
	int do_histogram, do_positivity, do_local_variation ;
	int do_blurring, do_bg_fitting, do_normalize_prtf ;
} ;

// algorithm.c
float DM_algorithm(struct algorithm_data*) ;
float HIO_algorithm(struct algorithm_data*) ;
float RAAR_algorithm(struct algorithm_data*) ;
float mod_DM_algorithm(struct algorithm_data*) ;
float ER_algorithm(struct algorithm_data*) ;

void make_recon_folders(struct algorithm_data*) ;
float run_iteration(struct algorithm_data*, int) ;
void save_current(struct algorithm_data*, int, struct timeval, struct timeval, float) ;
void save_output(struct algorithm_data*) ;
int parse_algorithm_strings(struct algorithm_data*, char*, char*) ;
void algorithm_allocate_memory(struct algorithm_data*) ;
void calc_prtf(struct algorithm_data*, int) ;
