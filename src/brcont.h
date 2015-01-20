#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <omp.h>

long size, vol, hklvol, num_supp ;
long hsize, ksize, lsize, hoffset, koffset, loffset ;
float *obs_mag, *exp_intens, *p1[3], *p2[3], *r1[3] ;
float *hkl_mag, *exp_hkl ;
long *support ;
fftwf_complex *fdensity, *rdensity ;
fftwf_complex *fhkl, *rhkl ;
fftwf_plan forward_cont, inverse_cont, forward_hkl, inverse_hkl ;

// diffmap.c
double diffmap(float**) ;

// setup.c
int setup() ;

// utils.c
void init_model(float**) ;
void average_model(float*, float*) ;
void gen_prtf(float*) ;
