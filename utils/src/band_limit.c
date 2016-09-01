#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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

int main(int argc, char *argv[]) {
	long i, size, vol, num_autocorr, num_supp ;
	float *intens ;
	uint8_t *support ;
	fftwf_complex *rdensity, *fdensity ;
	fftwf_plan forward, inverse ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <intens_fname> <size> <support_fname>\n", argv[0]) ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	vol = size*size*size ;
	
	// Read data
	intens = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	for (i = 0 ; i < vol ; ++i)
	if (intens[i] < -500.f)
		intens[i] = 0.f ;
	
	// Read support
	support = malloc(vol * sizeof(uint8_t)) ;
	fp = fopen(argv[3], "rb") ;
	fread(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	num_supp = 0 ;
	for (i = 0 ; i < vol ; ++i)
	if (support[i] == 1)
		num_supp++ ;
	
	// Initialize FFTW
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(omp_get_max_threads()) ;
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	
	sprintf(fname, "data/wisdom_%ld_%d", size, omp_get_max_threads()) ;
	fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_FORWARD, FFTW_MEASURE) ;
	}
	
	// Calculate autocorrelation of support
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = support[i] ;
	fftwf_execute(forward) ;
	for (i = 0 ; i < vol ; ++i)
		fdensity[i] = powf(cabsf(fdensity[i]), 2.f) / vol ;
	fftwf_execute(inverse) ;
	
	num_autocorr = 0 ;
	for (i = 0 ; i < vol ; ++i)
	if (crealf(rdensity[i]) > 0.5) {
		support[i] = 1 ;
		num_autocorr++ ;
	}
	else {
		support[i] = 0 ;
	}
	fprintf(stderr, "Number of voxels in support autocorrelation = %ld (%.2f)\n", num_autocorr, (double)num_autocorr/num_supp) ;
	
	// Constrain Patterson function of data by autocorrelation of support
	for (i = 0 ; i < vol ; ++i)
		fdensity[i] = intens[i] ;
	fftwf_execute(inverse) ;
	
	fp = fopen("data/rdensity.cpx", "wb") ;
	fwrite(rdensity, sizeof(fftwf_complex), vol, fp) ;
	fclose(fp) ;
	
	for (i = 0 ; i < vol ; ++i) {
		if (support[i] == 0)
			rdensity[i] = 0.f ;
		else
			rdensity[i] /= vol ;
	}
	fftwf_execute(forward) ;
	
	for (i = 0 ; i < vol ; ++i)
		intens[i] = crealf(fdensity[i]) ;
	
	sprintf(fname, "%s-smoothed.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving band-limited intensity to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(support) ;
	free(intens) ;
	free(rdensity) ;
	free(fdensity) ;
	
	return 0 ;
}
