#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "../../src/fft.h"
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long i, size, vol, num_autocorr, num_supp ;
	float *intens ;
	uint8_t *support ;
	FILE *fp ;
	char fname[999] ;
	struct fft_data fft ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <intens_fname> <support_fname>\n", argv[0]) ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
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
	fp = fopen(argv[2], "rb") ;
	fread(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	num_supp = 0 ;
	for (i = 0 ; i < vol ; ++i)
	if (support[i] == 1)
		num_supp++ ;
	
	// Initialize FFTW
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Calculate autocorrelation of support
	for (i = 0 ; i < vol ; ++i)
		fft.rdensity[i] = support[i] ;
	fft_forward(&fft) ;
	for (i = 0 ; i < vol ; ++i)
		fft.fdensity[i] = powf(cabsf(fft.fdensity[i]), 2.f) / vol ;
	fft_inverse(&fft) ;
	
	num_autocorr = 0 ;
	for (i = 0 ; i < vol ; ++i)
	if (crealf(fft.rdensity[i]) > 0.5) {
		support[i] = 1 ;
		num_autocorr++ ;
	}
	else {
		support[i] = 0 ;
	}
	fprintf(stderr, "Number of voxels in support autocorrelation = %ld (%.2f)\n", num_autocorr, (double)num_autocorr/num_supp) ;
	
	// Constrain Patterson function of data by autocorrelation of support
	for (i = 0 ; i < vol ; ++i)
		fft.fdensity[i] = intens[i] ;
	fft_inverse(&fft) ;
	
	fp = fopen("data/rdensity.cpx", "wb") ;
	fwrite(fft.rdensity, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	for (i = 0 ; i < vol ; ++i) {
		if (support[i] == 0)
			fft.rdensity[i] = 0.f ;
		else
			fft.rdensity[i] /= vol ;
	}
	fft_forward(&fft) ;
	
	for (i = 0 ; i < vol ; ++i)
		intens[i] = crealf(fft.fdensity[i]) ;
	
	sprintf(fname, "%s-smoothed.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving band-limited intensity to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(support) ;
	free(intens) ;
	fft_free(&fft) ;
	
	return 0 ;
}
