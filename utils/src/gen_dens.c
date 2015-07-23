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
	char fname[999] ;
	fftwf_plan inverse ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <cpx_fmodel> <size>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	c = size / 2 ;
	vol = size*size*size ;
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	temp = malloc(vol * sizeof(float)) ;
	
	sprintf(fname, "data/wisdom_%ld_16", size) ;
	fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Measuring plans\n") ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
		
		fp = fopen("data/wisdom_501_16", "wb") ;
		fftwf_export_wisdom_to_file(fp) ;
		fclose(fp) ;
	
		fprintf(stderr, "Created plans\n") ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}
	
	fp = fopen(argv[1], "rb") ;
	fread(rdensity, sizeof(fftwf_complex), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		fdensity[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
		  = rdensity[x*size*size + y*size + z] ;
	
	fftwf_execute(inverse) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
//		temp[((z+c)%size)*size*size + ((y+c)%size)*size + ((x+c)%size)] // Reversed axes
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

