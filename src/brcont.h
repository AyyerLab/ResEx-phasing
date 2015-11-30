//#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

#define IGNORE 2
#define BRAGG 1
#define NON_BRAGG 0

int iter ;
long size, vol, num_supp ;
float *iterate, *obs_mag, *p1, *p2, *r1 ;
float *coh_mag, *incoh_mag, *intens ;
long *support ;
uint8_t *mask ;
fftwf_complex *fdensity, *rdensity ;
float complex *bragg_calc ;
fftwf_plan forward, inverse ;

// diffmap.c
double diffmap(float*) ;
double error_red(float*) ;

// setup.c
int setup() ;
int setup_gen() ;

// utils.c
void init_model(float*) ;
void average_model(float*, float*) ;
void gen_prtf(float*) ;
void symmetrize_incoherent(fftwf_complex*, float*) ;
void symmetrize_coherent(fftwf_complex*, float*) ;
void blur_intens(float*, float*) ;
void apply_shrinkwrap(float*, float, float) ;
