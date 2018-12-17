#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

struct locval_pair {
	long loc ;
	float val ;
} ;

struct input_data {
	long size, num_supp, support_bounds[6] ;
	float *obs_mag ;
	uint8_t *support ;
	float complex *bragg_calc ;
	
	// do_histogram
	struct locval_pair *supp_locval ;
	float *inverse_cdf ;
	
	// do_local_variation
	struct locval_pair *local_variation ;
} ;

void input_init(struct input_data*, long) ;
void input_free(struct input_data*) ;

int input_parse_intens(struct input_data*, char*, float, int) ;
int input_parse_bragg(struct input_data*, char*, float) ;
int input_parse_support(struct input_data*, char*) ;
void input_init_iterate(struct input_data*, char*, char*, float*, int, int) ;
int input_read_histogram(struct input_data*, char*) ;

void input_match_histogram(struct input_data*, float*, float*) ;
void input_update_support(struct input_data*, float*, long) ;
void match_bragg(struct input_data*, float complex*, float) ;

