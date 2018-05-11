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
long size, vol, support_bounds[6] ;
float *obs_mag, *exp_mag ;
uint8_t *support ;
fftwf_complex *fdensity, *rdensity ;
float complex *bragg_calc ;
fftwf_plan forward_plan, inverse_plan ;
char output_prefix[999], point_group[999] ;

// Algorithm parameters
int num_iter, num_avg_iter ;
char (*algorithms)[8], (*avg_algorithms)[8] ;
int *intrad ;
double *radavg, *radcount, *obs_radavg ;
float algorithm_beta ;
float *algorithm_iterate, *algorithm_p1 ;
float *algorithm_p2, *algorithm_r1, *algorithm_r2 ;
int do_histogram, do_positivity, do_local_variation, do_normalize_prtf ;

// Radial background fitting
int do_bg_fitting ;

// Histogram matching
long num_supp, *supp_loc, *supp_index ;
float *supp_val, *inverse_cdf ;

// Support update
float *shrinkwrap_kernel ;
float *local_variation ;
long *voxel_pos ;

// Intensity blurring
int do_blurring, num_rot ;
double *quat ;

// Functions
// diffmap.c
float DM_algorithm(float*) ;
float HIO_algorithm(float*) ;
float RAAR_algorithm(float*) ;
float mod_DM_algorithm(float*) ;
float ER_algorithm(float*) ;
float modified_hio(float*) ;

// setup.c
int setup() ;
int setup_gen() ;

// utils.c
void init_model(float*, int, int) ;
void average_model(float*, float*) ;
void calc_prtf(float*, float*, int) ;
void symmetrize_incoherent(fftwf_complex*, float*, float*) ;
void blur_intens(float*, float*) ;
void apply_shrinkwrap(float*, float, float) ;
void dump_slices(float*, char*, int) ;
void dump_support_slices(uint8_t*, char*) ;
void init_radavg() ;
void radial_average(float*, float*) ;
void match_histogram(float*, float*) ;
float positive_mode(float*) ;
void variation_support(float*, uint8_t*, long) ;
void match_bragg(fftwf_complex*, float) ;
