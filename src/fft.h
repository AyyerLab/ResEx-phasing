#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/stat.h>
#include <libgen.h>

struct fft_data {
	long size ;
	char *wisdom_fname ;
	fftwf_complex *fdensity, *rdensity ;
	fftwf_plan forward_plan, inverse_plan ;
} ;

void fft_init(struct fft_data*, long, int) ;
void fft_create_plans(struct fft_data*) ;
void fft_gaussian_blur(struct fft_data*, float*, float) ;
void fft_apply_shrinkwrap(struct fft_data*, float*, float, float, uint8_t*, char*) ;
void fft_forward(struct fft_data*) ;
void fft_inverse(struct fft_data*) ;
void fft_free(struct fft_data*) ;
