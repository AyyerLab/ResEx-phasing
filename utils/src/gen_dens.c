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
		return 1 ;
	}
	size = atoi(argv[2]) ;
	c = size / 2 ;
	vol = size*size*size ;
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(4) ;
	
	fp = fopen("data/wisdom_501_4", "r") ;
	fftwf_import_wisdom_from_file(fp) ;
	fclose(fp) ;
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	temp = malloc(vol * sizeof(float)) ;
	
	inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(fdensity, sizeof(fftwf_complex), vol, fp) ;
	fclose(fp) ;
	
	fftwf_execute(inverse) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[((z+c)%size)*size*size + ((y+c)%size)*size + ((x+c)%size)]
		  = cabsf(rdensity[x*size*size + y*size + z]) / vol ;
	
	sprintf(fname, "%s_dens.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	fftwf_free(fdensity) ;
	fftwf_free(rdensity) ;
	fftwf_destroy_plan(inverse) ;
	
	return 0 ;
}

