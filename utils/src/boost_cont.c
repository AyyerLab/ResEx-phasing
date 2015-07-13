#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol ;
	float r ;
	float *model1 ;
	fftwf_complex *rdensity, *fdensity ;
	fftwf_plan forward, inverse ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <model1>\n", argv[0]) ;
		return 1 ;
	}
	size = 501 ;
	vol = size*size*size ;
	c = size / 2 ;
	
	// Allocate memory
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	model1 = malloc(vol * sizeof(float)) ;
	
	// Parse fftwf plan
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	if (size == 501) {
		fp = fopen("data/wisdom_501_16", "rb") ;
		if (fp == NULL) {
			fprintf(stderr, "Measuring plans\n") ;
			forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
			inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
			
			fp = fopen("data/wisdom_501_16", "wb") ;
			fftwf_export_wisdom_to_file(fp) ;
			fclose(fp) ;
			
			fprintf(stderr, "Created plans\n") ;
		}
		else {
			fftwf_import_wisdom_from_file(fp) ;
			fclose(fp) ;
			
			forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
			inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
		}
	}
	else 
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	
	// Parse model
	fp = fopen(argv[1], "rb") ;
	fread(model1, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate Fourier transform
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = model1[x] ;
	
	fftwf_execute(forward) ;
	
	for (x = 0 ; x < vol ; ++x)
		model1[x] = fdensity[x] ;
	
	fprintf(stderr, "Parsed model1\n") ;
	
	// Boost continuous transform
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		r = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		
		if (r > 0.65 * c)
			fdensity[x*size*size + y*size + z] *= 2. / vol ;
	}
	
	// Inverse Fourier transform
	fftwf_execute(inverse) ;
	
	for (x = 0 ; x < vol ; ++x)
		model1[x] = cabsf(rdensity[x]) ;
	
	// Write to file
	sprintf(fname, "%s-boost.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(model1, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model1) ;
	fftwf_free(rdensity) ;
	fftwf_free(fdensity) ;
	
	return 0 ;
}
