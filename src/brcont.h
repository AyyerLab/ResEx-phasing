//#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

int iter ;
long size, vol ;
float * restrict iterate, * restrict obs_mag, * restrict exp_mag, * restrict p1 ;
float *restrict p2, *restrict r1, *restrict r2 ;
uint8_t * restrict support ;
fftwf_complex * restrict fdensity, * restrict rdensity ;
float complex * restrict bragg_calc ;
fftwf_plan forward, inverse ;
char output_prefix[999], point_group[999] ;

// Histogram matching
long num_supp, *supp_loc, *supp_index ;
float *supp_val, *inverse_cdf ;

// diffmap.c
double diffmap(float*) ;
double error_red(float*) ;
double modified_hio(float*) ;

// setup.c
int setup() ;
int setup_gen() ;

// utils.c
void init_model(float*) ;
void average_model(float*, float*) ;
void gen_prtf(float*) ;
void symmetrize_incoherent(fftwf_complex*, float*) ;
void blur_intens(float*, float*) ;
void apply_shrinkwrap(float*, float, float) ;
void dump_slices(float*, char*) ;
void match_histogram(float*, float*) ;
