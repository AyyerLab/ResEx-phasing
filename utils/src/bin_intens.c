#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long x, y, z, s, size, red_size, factor ;
	float *model, *row ;
	int *counts ;
	char fname[500] ;
	FILE *fp  ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <model_fname> <shrink_factor>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	factor = strtol(argv[2], NULL, 10) ;
	red_size = (size / (2 * factor)) * 2 + 1 ;
	fprintf(stderr, "Reduced size = %ld\n", red_size) ;
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s_%ld.raw", remove_ext(argv[1]), red_size) ;
	
	row = malloc(size * sizeof(float)) ;
	model = calloc(red_size*red_size*red_size, sizeof(float)) ;
	counts = calloc(red_size*red_size*red_size, sizeof(int)) ;
	
	fp = fopen(argv[1], "rb") ;
	for (x = 0 ; x < size ; ++x) {
		if (x / factor >= red_size)
			continue ;
		
		for (y = 0 ; y < size ; ++y) {
			fread(row, sizeof(float), size, fp) ;
			
			if (y / factor >= red_size)
				continue ;
			
			for (z = 0 ; z < red_size ; ++z)
			for (s = 0 ; s < factor ; ++s)
			if (z*factor + s < size && row[z*factor + s] > 0.) {
				model[(x/factor)*red_size*red_size + (y/factor)*red_size + z] += row[z*factor + s] ;
				counts[(x/factor)*red_size*red_size + (y/factor)*red_size + z]++ ;
			}
		}
		if (x % 100 == 0)
			fprintf(stderr, "Finished x = %ld\n", x) ;
	}
	fclose(fp) ;
	
	for (x = 0 ; x < red_size*red_size*red_size ; ++x) {
		if (counts[x] > 0)
			model[x] /= counts[x] ;
		else
			model[x] = -1. ;
	}
	
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), red_size*red_size*red_size, fp) ;
	fclose(fp) ;
	
	float *temp = malloc(red_size*red_size*red_size*sizeof(float)) ;
	s = red_size ;
	long c = red_size / 2 ;
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z)
		temp[x*s*s + y*s + z] = 0.25 * (model[x*s*s + y*s + z] + 
		                                model[x*s*s + y*s + (2*c-z)] +
		                                model[x*s*s + (2*c-y)*s + z] +
	                                    model[x*s*s + (2*c-y)*s + (2*c-z)]) ;
	
	sprintf(fname, "%s_sym.raw", remove_ext(fname)) ;
	fp = fopen(fname, "wb") ;
	fwrite(temp, sizeof(float), red_size*red_size*red_size, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	free(counts) ;
	free(model) ;
	free(row) ;
	
	return 0 ;
}
