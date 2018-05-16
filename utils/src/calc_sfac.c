#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <omp.h>
#include "../../src/utils.h"

long size ;

long calc_bin(long x, long y, long z) {
	return (long) sqrt(x*x + y*y + z*z) ;
}

void calc_radavg(float *model, float *avg) {
	long x, y, z, c = size / 2, bin ;
	long *count = calloc(size, sizeof(long)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		bin = calc_bin(x-c, y-c, z-c) ;
		count[bin]++ ;
		avg[bin] += model[x*size*size + y*size + z] ;
	}
	
	for (x = 0 ; x < size ; ++x)
	if (count[x] > 0)
		avg[bin] /= count[bin] ;
	
	free(count) ;
}

float calc_sfac(float a, int n) {
	float numr = gsl_sf_gamma(n/2.+0.25)*gsl_sf_hyperg_1F1(n/2.+0.25, 0.5, a*a) + gsl_sf_gamma(n/2.+0.75)*gsl_sf_hyperg_1F1(n/2.+0.75, 1.5, a*a) ;
	float denr = gsl_sf_gamma(n/2.)*gsl_sf_hyperg_1F1(n/2., 0.5, a*a) + gsl_sf_gamma(n/2.+0.5)*gsl_sf_hyperg_1F1(n/2.+0.5, 1.5, a*a) ;
	
	return powf(2.f, 0.25f) * numr / denr ;
}

int main(int argc, char *argv[]) {
	long x, y, z, center, vol, i, bin ;
	float *intens, *sigma, *alpha, *radavg ;
	int num_twin = 4 ; // TODO Make general
	char fname[999] ;
	FILE *fp ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <intens_fname> <sigma_fname>\n", argv[0]) ;
		fprintf(stderr, "Optional: <output_fname>\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-sfac.raw", remove_ext(argv[1])) ;
	vol = size*size*size ;
	center = size / 2 ;
	
	intens = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed intens: %s\n", argv[1]) ;
	
	sigma = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[2], "rb") ;
	fread(sigma, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed sigma: %s\n", argv[2]) ;
	
	radavg = calloc(size, sizeof(float)) ;
	calc_radavg(intens, radavg) ;
	fprintf(stderr, "Calculated radial average\n") ;
	
	alpha = calloc(vol, sizeof(float)) ;
	#pragma omp parallel default(shared) private(x,y,z,i,bin)
	{
		int rank = omp_get_thread_num() ;
		
		#pragma omp for schedule(static,1)
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			i = x*size*size + y*size + z ;
			if (intens[i] > -5000.f) {
				bin = calc_bin(x-center, y-center, z-center) ;
				alpha[i] = (intens[i]*radavg[bin] - num_twin*sigma[i]*sigma[i]) / (sqrtf(2.f) * radavg[bin] * sigma[i]) ;
			}
			if (rank == 0 && y == size-1 && z == size-1)
				fprintf(stderr, "\rCalculated alpha for x = %ld", x) ;
		}
	}
	fprintf(stderr, "\n") ;
	
/*	fp = fopen("data/ps2-qdark-alpha.raw", "wb") ;
	fwrite(alpha, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	return 1 ;
*/	
	#pragma omp parallel default(shared) private(i)
	{
		int rank = omp_get_thread_num() ;
		
		#pragma omp for schedule(static,1)
		for (i = 0 ; i < vol ; ++i) {
			if (intens[i] < -5000.f)
				intens[i] = -1.f ;
			else if (sigma[i] == 0.f)
				intens[i] = -1.f ;
			else
				intens[i] = sqrtf(sigma[i]) * calc_sfac(alpha[i], num_twin) ;
			
			if (rank == 0 && i % (size*size) == size*size - 1)
				fprintf(stderr, "\rCalculated sfac for x = %ld", i/(size*size)) ;
		}
	}
	fprintf(stderr, "\n") ;
	
	fprintf(stderr, "Writing output to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(sigma) ;
	free(intens) ;
	free(alpha) ;
	free(radavg) ;
	
	return 0 ;
}
