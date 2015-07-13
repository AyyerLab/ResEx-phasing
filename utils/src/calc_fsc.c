#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, bin, num_bins, vox ;
	double d_max, r, binsize ;
	float *temp, *norm1, *norm2 ;
	float complex *model1, *model2, *fsc ;
	fftwf_complex *rdensity, *fdensity ;
	fftwf_plan forward ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <model1> <model2>\n", argv[0]) ;
		return 1 ;
	}
	size = 501 ;
	vol = size*size*size ;
	c = size / 2 ;
	
	// Allocate memory
	temp = malloc(vol * sizeof(float)) ;
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	model1 = malloc(vol * sizeof(float complex)) ;
	model2 = malloc(vol * sizeof(float complex)) ;
	
	num_bins = 50 ;
	d_max = 2200./750. ;
	binsize = (double) c / num_bins + 1 ;
	
	fsc = calloc(num_bins, sizeof(float complex)) ;
	norm1 = calloc(num_bins, sizeof(float)) ;
	norm2 = calloc(num_bins, sizeof(float)) ;
	
	// Parse fftwf plan
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	if (size == 501) {
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
	}
	else 
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;
	
	// Parse first model
	fp = fopen(argv[1], "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate Fourier transform
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = temp[x] ;
	
	fftwf_execute(forward) ;
	
	for (x = 0 ; x < vol ; ++x)
		model1[x] = fdensity[x] ;
	
	fprintf(stderr, "Parsed model1\n") ;
	
	// Parse second model
	fp = fopen(argv[2], "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate Fourier transform
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = temp[x] ;
	
	fftwf_execute(forward) ;
	
	for (x = 0 ; x < vol ; ++x)
		model2[x] = fdensity[x] ;
	
	fprintf(stderr, "Parsed model2\n") ;
	
	// Calculate FSC
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		r = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		bin = (int) (r / binsize) ;
		if (bin > num_bins)
			continue ;
		
		vox = x*size*size + y*size + z ;
		fsc[bin] += model1[vox] * conjf(model2[vox]) ;
		norm1[bin] += powf(cabsf(model1[vox]), 2.f) ;
		norm2[bin] += powf(cabsf(model2[vox]), 2.f) ;
	}
	
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (norm1[bin] * norm2[bin] > 0.)
		fsc[bin] /= sqrtf(norm1[bin] * norm2[bin]) ;
	
	// Write to file
	sprintf(fname, "%s", extract_fname(argv[1])) ;
	strtok(fname, "_.") ;
	int num1 = atoi(strtok(NULL, "_.")) ;
	sprintf(fname, "%s", extract_fname(argv[2])) ;
	strtok(fname, "_.") ;
	int num2 = atoi(strtok(NULL, "_.")) ;
	sprintf(fname, "logs/fsc-%d-%d.dat", num1, num2) ;
	fprintf(stderr, "Writing to %s\n", fname) ;
	
	fp = fopen(fname, "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		fprintf(fp, "%.3f\t%.3f\t%.6f\n", (bin+1)/d_max/num_bins, num_bins*d_max/(bin+1), cabsf(fsc[bin])) ;
	fclose(fp) ;
	/*
	// Free memory
	free(temp) ;
	free(model1) ;
	free(model2) ;
	fftwf_free(rdensity) ;
	fftwf_free(fdensity) ;
	free(fsc) ;
	free(norm1) ;
	free(norm2) ;
	*/
	return 0 ;
}
