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

int iter ;
long size, vol, hklvol, num_supp ;
long hsize, ksize, lsize, hoffset, koffset, loffset ;
float *iterate, *obs_mag, *exp_mag, *p1, *p2, *r1 ;
long *support ;
int num_rot ;
double *quat ;
fftwf_complex *fdensity, *rdensity ;
fftwf_complex *fhkl, *rhkl, *hkl_calc ;
fftwf_plan forward_cont, inverse_cont, forward_hkl, inverse_hkl ;

// diffmap.c
double diffmap(float*) ;

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
