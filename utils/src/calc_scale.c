#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <math.h>
#include <locale.h>
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	long x, y, z, vox ;
	long size, c, vol, num_vox = 0 ;
	double dot = 0., normsq = 0., rsq, scale, rmin, rmax ;
	struct ccp4_map obs_map = {0}, model_map = {0} ;
	setlocale(LC_ALL, "C") ; // For commas in large integers
	
	if (argc < 5) {
		fprintf(stderr, "Calc Scale: Calculates scale factor between two intensities\n") ;
		fprintf(stderr, "-----------------------------------------------------------\n") ;
		fprintf(stderr, "Compares two intensity volumes in a given radius range\n") ;
		fprintf(stderr, "\nUsage: %s <sym_model_fname> <merge_fname> <rmin> <rmax>\n", argv[0]) ;
		return 1 ;
	}
	rmin = atof(argv[3]) ;
	rmax = atof(argv[4]) ;
	
	// Parse sym_model and obs_map
	if (parse_map(argv[1], &model_map))
		return 1 ;
	if (model_map.f32_data == NULL) {
		fprintf(stderr, "sym_model needs float data\n") ;
		return 1 ;
	}
	if (parse_map(argv[2], &obs_map))
		return 1 ;
	if (obs_map.f32_data == NULL) {
		fprintf(stderr, "merge %s needs float data\n", argv[2]) ;
		return 1 ;
	}
	
	size = model_map.header.nx ;
	c = size / 2 ;
	vol = size*size*size ;
	
	// Calculate magnitude
	for (x = 0 ; x < vol ; ++x) {
		model_map.f32_data[x] = sqrtf(model_map.f32_data[x]) ;
		if (obs_map.f32_data[x] > 0.f)
			obs_map.f32_data[x] = sqrtf(obs_map.f32_data[x]) ;
	}
	
	// Calculate scale factor
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rsq = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c) ;
		vox = x*size*size + y*size + z ;
		
		if (rsq > rmin*rmin && rsq < rmax*rmax && obs_map.f32_data[vox] > 0.) {
			//dot += model_map.f32_data[vox] ;
			//normsq += obs_map.f32_data[vox] ;
			dot += obs_map.f32_data[vox] * model_map.f32_data[vox] ;
			normsq += pow(obs_map.f32_data[vox], 2.) ;
			num_vox++ ;
		}
	}
	
	scale = dot / normsq ;
	
	printf("Scale factor for magnitude: %f\n", scale) ;
	printf("%'ld voxels contributed to this calculation\n", num_vox) ;
	
	// Free memory
	free_map(&model_map) ;
	free_map(&obs_map) ;
	
	return 0 ;
}
