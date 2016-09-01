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
	long x, y, z, size, c, vol ;
	fftwf_complex *fdensity, *rdensity ;
	float *temp, r_max = -1, dist ;
	FILE *fp ;
	char fname[999] ;
	fftwf_plan inverse ;
	
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
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(32) ;
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	temp = malloc(vol * sizeof(float)) ;
	
	sprintf(fname, "data/wisdom_%ld_32", size) ;
	fp = fopen(fname, "rb") ;
	if (fp == NULL)
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_ESTIMATE) ;
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}
	
	// Read complex Fourier amplitudes
	fp = fopen(argv[1], "rb") ;
	fread(rdensity, sizeof(fftwf_complex), vol, fp) ;
	fclose(fp) ;
	
	// Shift coordinates such that origin is at the corner
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		if (r_max < 0.)
			fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
			  = rdensity[x*size*size + y*size + z] ;
		else {
			dist = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
			if (dist > r_max)
				fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = 0. ;
			else
				fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
				  = rdensity[x*size*size + y*size + z] ;
		}
	}
	
	// Do inverse Fourier transform
	fftwf_execute(inverse) ;
	
	// Shift real space density such that it is in the center of the cube
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
		  = crealf(rdensity[x*size*size + y*size + z]) / vol ;
	
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-dens.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	fftwf_free(fdensity) ;
	fftwf_free(rdensity) ;
	fftwf_destroy_plan(inverse) ;
	
	return 0 ;
}

