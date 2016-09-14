#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
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
	float *intens, *radavg, *radsqavg, *sigma, *kl_div ;
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
	radavg = calloc(size, sizeof(float)) ;
	radsqavg = calloc(size, sizeof(float)) ;
	sigma = malloc(size * sizeof(float)) ;
	kl_div = malloc(size * sizeof(float)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed intens: %s\n", argv[1]) ;
	
	for (i = 0 ; i < vol ; ++i)
	if (intens[i] < -500.)
		intens[i] = 0. ;
	
	rmax = calc_radavg(intens, radavg) ;
	fprintf(stderr, "Calculated radial average\n") ;
	
	for (i = 0 ; i < vol ; ++i)
		intens[i] = intens[i]*intens[i] ;
	
	calc_radavg(intens, radsqavg) ;
	fprintf(stderr, "Calculated squared radial average\n") ;
	
	int num_sym = 4 ;
	for (i = 0 ; i < size ; ++i) {
		sigma[i] = radavg[i] / sqrtf(radsqavg[i] - radavg[i]*radavg[i]) ;
		if (isnan(sigma[i]))
			sigma[i] = 0.f ;
		
		kl_div[i] = num_sym * log(num_sym / radavg[i])
		          + log(sqrt(2.*M_PI) * sigma[i] / gsl_sf_fact(num_sym-1))
		          + (num_sym-1)*(gsl_sf_psi_int(num_sym) - log(num_sym/radavg[i]))
				  - num_sym
				  + radavg[i]*radavg[i] / (2. * sigma[i]*sigma[i]) ;
		if (isnan(kl_div[i]))
			kl_div[i] = 0.f ;
	}
	
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-snr.dat", remove_ext(argv[1])) ;
	fprintf(stderr, "Writing output to %s\n", fname) ;
	fp = fopen(fname, "w") ;
	for (i = 0 ; i < rmax ; ++i)
		fprintf(fp, "%.4ld\t%.6e\t%.6e\t%.6e\n", i, radavg[i], sigma[i], kl_div[i]) ;
	fclose(fp) ;
	
	free(intens) ;
	free(radavg) ;
	free(radsqavg) ;
	
	return 0 ;
}
