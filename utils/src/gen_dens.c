#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "../../src/fft.h"

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol ;
	float *temp, r_max = -1, dist ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <cpx_fmodel> <size>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "Second option: <r_max> cutoff radius\n") ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	c = size / 2 ;
	vol = size*size*size ;
	
	if (argc > 4) {
		r_max = atof(argv[4]) ;
		fprintf(stderr, "Truncated to radius = %f\n", r_max) ;
	}
	
	struct fft_data fft ;
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	temp = malloc(vol * sizeof(float)) ;
	
	// Read complex Fourier amplitudes
	fp = fopen(argv[1], "rb") ;
	fread(fft.rdensity, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Shift coordinates such that origin is at the corner
	// Also truncate to radius if specified
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		if (r_max < 0.)
			fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
			  = fft.rdensity[x*size*size + y*size + z] ;
		else {
			dist = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
			if (dist > r_max)
				fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = 0. ;
			else
				fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
				  = fft.rdensity[x*size*size + y*size + z] ;
		}
	}
	
	// Do inverse Fourier transform
	fft_inverse(&fft) ;
	
	// Shift real space density such that it is in the center of the cube
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
		  = crealf(fft.rdensity[x*size*size + y*size + z]) / vol ;
	
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-dens.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	fft_free(&fft) ;
	
	return 0 ;
}

