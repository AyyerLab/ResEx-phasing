#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <omp.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, num_bins, vox, *numvox ;
	double d_min, r, binsize ;
	int *radbin ;
	double *norm1, *norm2, *cc, *radavg1, *radavg2 ;
	float *intens1, *intens2 ;
	FILE *fp ;
	char fname[999] ;

	if (argc < 4) {
		fprintf(stderr, "Format: %s <intens1> <intens2> <res_at_edge>\n", argv[0]) ;
		fprintf(stderr, "Optional: <cc_fname>\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	d_min = atof(argv[3]) ;
	if (argc > 4) {
		strcpy(fname, argv[4]) ;
	}
	else {
		sprintf(fname, "%s", extract_fname(argv[1])) ;
		strtok(fname, "_.") ;
		int num1 = atoi(strtok(NULL, "_.")) ;
		sprintf(fname, "%s", extract_fname(argv[2])) ;
		strtok(fname, "_.") ;
		int num2 = atoi(strtok(NULL, "_.")) ;
		sprintf(fname, "cc-%d-%d.dat", num1, num2) ;
	}
	
	vol = size*size*size ;
	c = size / 2 ;
	binsize = 1. ;
	num_bins = ceil(c / binsize) + 1 ;
	fprintf(stderr, "Model size = %ld, num_bins = %ld\n", size, num_bins) ;

	// Allocate memory
	intens1 = malloc(vol * sizeof(float)) ;
	intens2 = malloc(vol * sizeof(float)) ;
	radbin = malloc(vol * sizeof(int)) ;
	
	cc = calloc(num_bins, sizeof(double)) ;
	norm1 = calloc(num_bins, sizeof(double)) ;
	norm2 = calloc(num_bins, sizeof(double)) ;
	numvox = calloc(num_bins, sizeof(long)) ;
	radavg1 = calloc(num_bins, sizeof(double)) ;
	radavg2 = calloc(num_bins, sizeof(double)) ;

	// Parse intensities
	fp = fopen(argv[1], "rb") ;
	fread(intens1, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fp = fopen(argv[2], "rb") ;
	fread(intens2, sizeof(float), vol, fp) ;
	fclose(fp) ;

	// Subtract radial averages
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		vox = x*size*size + y*size + z ;
		r = sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		radbin[vox] = (int) (r / binsize) ;
		if (radbin[vox] > num_bins - 1)
			continue ;
		
		radavg1[radbin[vox]] += intens1[vox] ;
		radavg2[radbin[vox]] += intens2[vox] ;
		numvox[radbin[vox]]++ ;
	}
	
	for (x = 0 ; x < num_bins ; ++x)
	if (numvox[x] > 0) {
		radavg1[x] /= numvox[x] ;
		radavg2[x] /= numvox[x] ;
	}
	
	for (vox = 0 ; vox < vol ; ++vox)
	if (radbin[vox] < num_bins) {
		intens1[vox] -= radavg1[radbin[vox]] ;
		intens2[vox] -= radavg2[radbin[vox]] ;
	}

	// Calculate CC between two intensities
	for (vox = 0 ; vox < vol ; ++vox) {
		if (radbin[vox] > num_bins - 1)
			continue ;
		
		cc[radbin[vox]] += intens1[vox] * intens2[vox] ;
		norm1[radbin[vox]] += pow(intens1[vox], 2.) ;
		norm2[radbin[vox]] += pow(intens2[vox], 2.) ;
	}
	
	for (x = 0 ; x < num_bins ; ++x)
	if (norm1[x] * norm2[x] > 0.)
		cc[x] /= sqrt(norm1[x] * norm2[x]) ;

	// Write to file
	fprintf(stderr, "Writing to %s\n", fname) ;
	fp = fopen(fname, "w") ;
	for (x = 0 ; x < num_bins ; ++x)
		fprintf(fp, "%6.3f  %6.3f  %8.6f  %.8ld\n", 
		        (x+1)/d_min/num_bins, 
		        num_bins*d_min/(x+1), 
		        cc[x], 
		        numvox[x]) ;
	fclose(fp) ;

	// Free memory
	free(intens1) ;
	free(intens2) ;
	free(radbin) ;
	free(cc) ;
	free(norm1) ;
	free(norm2) ;
	free(radavg1) ;
	free(radavg2) ;
	free(numvox) ;
	
	return 0 ;
}
