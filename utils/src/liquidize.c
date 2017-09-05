#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <omp.h>

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

long factorial(long i) {
	// Messy factorial function
	if (i < 2)
		return i ;
	return i*factorial(i - 1) ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, fac ;
	float vox_size, qscale, sigma, gamma ;
	double error, mean ;
	float *radius, *intens, *autocorr ;
	fftwf_complex *rdensity, *fdensity ;
	fftwf_plan forward, inverse ;
	FILE *fp ;
	char fname[1024] ;
	
	if (argc < 6) {
		fprintf(stderr, "Format: %s <cpx_fname> <size> <res_at_edge> <sigma> <gamma>\n", argv[0]) ;
		fprintf(stderr, "<res_at_edge>, <sigma> and <gamma> have consistent length units\n") ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	c = size / 2 ;
	vol = size*size*size ;
	vox_size = atof(argv[3]) / 2. ;
	qscale = 1. / (2. * vox_size * c) ;
	sigma = atof(argv[4]) ;
	gamma = atof(argv[5]) ;

	// Allocate memory and FFTW plans
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	autocorr = malloc(vol * sizeof(float)) ;
	intens = calloc(vol, sizeof(float)) ;
	radius = malloc(vol * sizeof(float)) ;
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(32) ;
	sprintf(fname, "data/wisdom_%ld_32", size) ;
	fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_ESTIMATE) ;
	}
	else {
		fprintf(stderr, "Reading wisdom from %s\n", fname) ;
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}

	// Read complex Fourier amplitudes
	fp = fopen(argv[1], "rb") ;
	fread(rdensity, sizeof(fftwf_complex), vol, fp) ;
	fclose(fp) ;

	// Calculate intensities and shift origin
	mean = 0. ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
		  = powf(cabsf(rdensity[x*size*size + y*size + z]), 2.f) ;
		radius[x*size*size + y*size + z] = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		mean += fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] ;
	}

	// Get ideal autocorrelation
	fftwf_execute(inverse) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		autocorr[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)] = crealf(rdensity[x*size*size + y*size + z]) / vol ;

	// Calculate modified intensities
	fprintf(stderr, "Starting convolutions\n") ;
	float *diff_array = malloc(vol * sizeof(float)) ;
	int break_flag = 0 ;
	#pragma omp parallel default(shared)
	{
		int n, omp_rank = omp_get_thread_num() ;
		long x, y, z, vox ;
		double diff, exponent ;
		for (n = 1 ; n <= 100 ; ++n) {
			if (break_flag)
				continue ;
			
			#pragma omp for schedule(static,1)
			for (x = 0 ; x < size ; ++x)
			for (y = 0 ; y < size ; ++y)
			for (z = 0 ; z < size ; ++z) {
				vox = x*size*size + y*size + z ;
				exponent = n * radius[vox] * vox_size / gamma ;
				rdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
				  = autocorr[vox] * expf(-exponent) ;
			}
			
			if (omp_rank == 0) {
				fftwf_execute(forward) ;
				fac = factorial(n) ;
				error = 0. ;
			}
			#pragma omp barrier
			
			#pragma omp for schedule(static,1) reduction(+:error)
			for (x = 0 ; x < size ; ++x)
			for (y = 0 ; y < size ; ++y)
			for (z = 0 ; z < size ; ++z) {
				vox = x*size*size + y*size + z ;
				exponent = powf(2. * M_PI * sigma * radius[vox] * qscale, 2.f) ;
				if (n < 10) {
					diff = cabs(fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]) *
						   pow(exponent, n) / fac * exp(-exponent) ;
				}
				else { // Stirling's approximation for n >= 10
					diff = log(cabs(fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)])) +
						   n*log(exponent) - (log(2.*M_PI*n)/2 + n*log(n) - n) - exponent ;
					diff = exp(diff) ;
				}
				if (radius[vox] > 0.)
					diff /= 1. - exp(-exponent) ;
				diff_array[vox] = diff ;
				error += diff ;
				intens[vox] += diff ;
			}
			
			if (omp_rank == 0) {
				sprintf(fname, "data/liq-%.2d.raw", n) ;
				fp = fopen(fname, "wb") ;
				fwrite(diff_array, sizeof(float), vol, fp) ;
				fclose(fp) ;
				
				fprintf(stderr, "Finished n = %.3d, rel_error = %e\n", n, error/mean) ;
				if (n > 2 && error/mean < 1.e-5)
					break_flag = 1 ;
			}
			#pragma omp barrier
		}
	}

	// Save output
	sprintf(fname, "%s-liq.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving output to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;

	// Free memory
	free(rdensity) ;
	free(fdensity) ;
	free(autocorr) ;
	free(intens) ;
	free(radius) ;

	return 0 ;
}
