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
	float *temp ;
	FILE *fp ;
	char fname[999], symfname[999] ;
	fftwf_plan forward ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <raw_model> <size>\n", argv[0]) ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	c = size / 2 ;
	vol = size*size*size ;
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	// Allocate memory
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	temp = malloc(vol * sizeof(float)) ;
	
	// Generate FFTW plans
	sprintf(fname, "data/wisdom_%ld_32", size) ;
	fp = fopen(fname, "rb") ;
	if (fp == NULL)
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
	}
	
	// Parse real-space density
	fp = fopen(argv[1], "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Translate array such that molecule is around the origin
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
//		rdensity[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)] = temp[x*size*size + y*size + z] ;
		rdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = temp[x*size*size + y*size + z] ;
	
	// Apply Fourier transform
	fftwf_execute(forward) ;
	
	// Translate array to put the origin at (c, c, c)
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		rdensity[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
		  = fdensity[x*size*size + y*size + z] ; 
	
	// Write complex array to file
	if (argc > 3) {
		strcpy(fname, argv[3]) ;
		sprintf(symfname, "%s-sym.raw", remove_ext(fname)) ;
	}
	else {
		sprintf(fname, "%s-fdens.cpx", remove_ext(argv[1])) ;
		sprintf(symfname, "%s-fdens-sym.raw", remove_ext(argv[1])) ;
	}
	fp = fopen(fname, "wb") ;
	fwrite(rdensity, sizeof(fftwf_complex), vol, fp) ;
	fclose(fp) ;
	
	// Symmetrize intensities
	memset(temp, 0, vol * sizeof(float)) ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[x*size*size + y*size + z] = 0.25 * (pow(cabsf(rdensity[x*size*size + y*size + z]), 2.) + 
		                                         pow(cabsf(rdensity[x*size*size + y*size + (2*c-z)]), 2.) +
		                                         pow(cabsf(rdensity[x*size*size + (2*c-y)*size + z]), 2.) +
		                                         pow(cabsf(rdensity[x*size*size + (2*c-y)*size + (2*c-z)]), 2.)) ;
	
	// Write symmetrized intensities to file
	fp = fopen(symfname, "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(temp) ;
	fftwf_free(fdensity) ;
	fftwf_free(rdensity) ;
	fftwf_destroy_plan(forward) ;
	
	return 0 ;
}

