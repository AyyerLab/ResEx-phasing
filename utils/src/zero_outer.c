#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol ;
	float rad, rmin, rmax, minmodel = FLT_MAX ;
	float *model ;
	FILE *fp ;
	char fname[500] ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <model_fname> <rmin> <rmax>\n", argv[0]) ;
		fprintf(stderr, "Optional: <subtract_min>. Default 0.\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	rmin = atof(argv[2]) ;
	rmax = atof(argv[3]) ;
	
	vol = size*size*size ;
	c = size / 2 ;
	
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Read model\n") ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rad = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		
		if (rad < rmin)
			model[x*size*size + y*size + z] = -1.f ;
		else if (rad > rmax + 3)
			model[x*size*size + y*size + z] = 0.f ;
		else if (rad > rmax - 3)
			model[x*size*size + y*size + z] = -1.f ;
		else if (model[x*size*size + y*size + z] < minmodel && model[x*size*size + y*size + z] > -1000.)
			minmodel = model[x*size*size + y*size + z] ;
	}
	fprintf(stderr, "Zeroed outer voxels and set inner voxels to be negative\n") ;
	fprintf(stderr, "Minimum value in annulus = %.3e\n", minmodel) ;
	
	if (argc > 4 && atoi(argv[4]) > 0) {
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			rad = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
			
			if (rad >= rmin && rad <= rmax - 3)
				model[x*size*size + y*size + z] -= minmodel ;
		}
		fprintf(stderr, "Subtracted minimum value in annulus\n") ;
		sprintf(fname, "%s-pos-zero.raw", remove_ext(argv[1])) ;
	}
	else {
		sprintf(fname, "%s-zero.raw", remove_ext(argv[1])) ;
	}
	fprintf(stderr, "%s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	
	return 0 ;
}
