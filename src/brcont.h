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
#include <sys/stat.h>

int iter ;
long size, vol ;
float * restrict obs_mag, * restrict exp_mag ;
uint8_t * restrict support ;
fftwf_complex * restrict fdensity, * restrict rdensity ;
float complex * restrict bragg_calc ;
fftwf_plan forward_plan, inverse_plan ;
char output_prefix[999], point_group[999] ;

// Algorithm parameters
char algorithm_name[999], avg_algorithm_name[999] ;
int *intrad ;
float *radavg, *radcount, algorithm_beta ;
float * restrict algorithm_iterate, * restrict algorithm_p1 ;
float * restrict algorithm_p2, * restrict algorithm_r1, * restrict algorithm_r2 ;
int do_histogram, do_positivity ;

// Histogram matching
long num_supp, *supp_loc, *supp_index ;
float *supp_val, *inverse_cdf ;

// Functions
// diffmap.c
double DM_algorithm(float*) ;
double HIO_algorithm(float*) ;
double RAAR_algorithm(float*) ;
double mod_DM_algorithm(float*) ;
double ER_algorithm(float*) ;
double modified_hio(float*) ;

// setup.c
int setup() ;
int setup_gen() ;

// utils.c
void init_model(float*, int) ;
void average_model(float*, float*) ;
void gen_prtf(float*) ;
void symmetrize_incoherent(fftwf_complex*, float*) ;
void blur_intens(float*, float*) ;
void apply_shrinkwrap(float*, float, float) ;
void dump_slices(float*, char*, int) ;
void init_radavg() ;
void radial_average(float*, float*) ;
void match_histogram(float*, float*) ;
float positive_mode(float*) ;
