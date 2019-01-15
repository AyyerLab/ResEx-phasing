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
		fprintf(stderr, "Bin Intens: Downsample given volume by integer factor\n") ;
		fprintf(stderr, "-----------------------------------------------------\n") ;
		fprintf(stderr, "Reduces size of 3D volume by combining nxnxn voxels\n") ;
		fprintf(stderr, "\nUsage: %s <model_fname> <shrink_factor>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "\nOutput: <model_fname>-<reduced_size>.raw (if <out_fname> not given)\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	factor = strtol(argv[2], NULL, 10) ;
	red_size = (size / (2 * factor)) * 2 + 1 ;
	fprintf(stderr, "Reduced size = %ld\n", red_size) ;
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-%ld.raw", remove_ext(argv[1]), red_size) ;
	
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
	
	fprintf(stderr, "Saving output to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), red_size*red_size*red_size, fp) ;
	fclose(fp) ;
	
	free(counts) ;
	free(model) ;
	free(row) ;
	
	return 0 ;
}
