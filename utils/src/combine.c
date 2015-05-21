#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	long x, y, z, i, size, c, vol ;
	fftwf_complex *fdensity, *rdensity ;
	float *recon, *temp ;
	FILE *fp ;
	char fname[999] ;
	fftwf_plan forward ;
	
	if (argc < 3) {
		fprintf(stderr, "Averages a set of reconstructions in the results folder\n") ;
		fprintf(stderr, "Format: %s <num1> <num2> ...\n", argv[0]) ;
		return 1 ;
	}
	size = 501 ;
	c = size / 2 ;
	vol = size*size*size ;
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	recon = malloc(vol * sizeof(float)) ;
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	temp = malloc(vol * sizeof(float)) ;
	
	// Parse FFTW wisdom
	fp = fopen("data/wisdom_501_16", "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Measuring plans\n") ;
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		
		fp = fopen("data/wisdom_501_16", "wb") ;
		fftwf_export_wisdom_to_file(fp) ;
		fclose(fp) ;
	
		fprintf(stderr, "Created plans\n") ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
	}
	
	// Calculate mean real-space density
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = 0.f ;
	
	for (i = 1 ; i < argc ; ++i) {
		sprintf(fname, "results/output_%.2d.raw", atoi(argv[i])) ;
		fp = fopen(fname, "rb") ;
		fread(recon, sizeof(float), vol, fp) ;
		fclose(fp) ;
		
		for (x = 0 ; x < vol ; ++x)
			rdensity[x] += recon[x] ;
	}
	
	if (argc > 1)
	for (x = 0 ; x < vol ; ++x) {
		rdensity[x] /= (argc - 1) ;
		temp[x] = crealf(rdensity[x]) ;
	}
	
	// Write mean real-space density to file
	fp = fopen("results/comb_output.raw", "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate mean Fourier intensity
	fftwf_execute(forward) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[((z+c)%size)*size*size + ((y+c)%size)*size + ((x+c)%size)]
		  = cabsf(fdensity[x*size*size + y*size + z]) ;
	
	fp = fopen("results/comb_foutput.raw", "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	free(recon) ;
	fftwf_free(fdensity) ;
	fftwf_free(rdensity) ;
	fftwf_destroy_plan(forward) ;
	
	return 0 ;
}

