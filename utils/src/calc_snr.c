#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <omp.h>

long size ;

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

double calc_bin(long x, long y, long z) {
	return floor(sqrt(x*x + y*y + z*z)) ;
}

long calc_radavg(float *model, float *avg) {
	long x, y, z, c = size / 2, bin, rmax = 0 ;
	long *count = calloc(size, sizeof(long)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		bin = (long) calc_bin(x-c, y-c, z-c) ;
		count[bin]++ ;
		avg[bin] += model[x*size*size + y*size + z] ;
	}
	
	for (x = 0 ; x < size ; ++x)
	if (count[x] > 0) {
		avg[x] /= count[x] ;
		rmax = x ;
	}
	
	free(count) ;
	
	return rmax ;
}

int main(int argc, char *argv[]) {
	long vol, i, rmax ;
	float *intens, *radavg, *radsqavg ;
	char fname[999] ;
	FILE *fp ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <intens_fname> <size>\n", argv[0]) ;
		fprintf(stderr, "Optional: <output_fname>\n") ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	vol = size*size*size ;
	
	intens = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed intens: %s\n", argv[1]) ;
	
	for (i = 0 ; i < vol ; ++i)
	if (intens[i] < -500.)
		intens[i] = 0. ;
	
	radavg = calloc(size, sizeof(float)) ;
	rmax = calc_radavg(intens, radavg) ;
	fprintf(stderr, "Calculated radial average\n") ;
	
	for (i = 0 ; i < vol ; ++i)
		intens[i] = intens[i]*intens[i] ;
	
	radsqavg = calloc(size, sizeof(float)) ;
	calc_radavg(intens, radsqavg) ;
	fprintf(stderr, "Calculated squared radial average\n") ;
	
	// The factor '5/4' in the following fomula is for num_twin = 4
	for (i = 0 ; i < size ; ++i) {
//		radavg[i] = radavg[i] / sqrtf(radsqavg[i] - 1.25 * radavg[i]*radavg[i]) ;
		radavg[i] = radavg[i] / sqrtf(radsqavg[i] - radavg[i]*radavg[i]) ;
		if (isnan(radavg[i]))
			radavg[i] = 0.f ;
	}
	
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-snr.dat", remove_ext(argv[1])) ;
	fprintf(stderr, "Writing output to %s\n", fname) ;
	fp = fopen(fname, "w") ;
	for (i = 0 ; i < rmax ; ++i)
		fprintf(fp, "%.4ld\t%.6e\n", i, radavg[i]) ;
	fclose(fp) ;
	
	free(intens) ;
	free(radavg) ;
	free(radsqavg) ;
	
	return 0 ;
}
