#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	long x, y, z, size, c ;
	float rad, rmin, rmax, minmodel = FLT_MAX ;
	char fname[1024] ;
	struct ccp4_map map = {0} ;
	
	if (argc < 4) {
		fprintf(stderr, "Zero Outer: Process outer and inner parts of intensity file\n") ;
		fprintf(stderr, "-----------------------------------------------------------\n") ;
		fprintf(stderr, "voxels outside <rmax> are zeroed and ones inside <rmin> are set to -1.\n") ;
		fprintf(stderr, "If <subtract_min> is non-zero, the minimum value in each radial shell is subtracted\n") ;
		fprintf(stderr, "\nUsage: %s <intens_fname> <rmin> <rmax>\n", argv[0]) ;
		fprintf(stderr, "\nOutput: <intens_fname>-zero.ccp4\n") ;
		return 1 ;
	}
	rmin = atof(argv[2]) ;
	rmax = atof(argv[3]) ;
	
	if (parse_map(argv[1], &map))
		return 1 ;
	if (map.f32_data == NULL) {
		fprintf(stderr, "Need float data in input intens\n") ;
		return 1 ;
	}
	if (map.header.nx != map.header.ny || map.header.nx != map.header.nz) {
		fprintf(stderr, "Need cubic volume of data (input shape: %d, %d, %d)\n", map.header.nx, map.header.ny, map.header.nz) ;
		return 1 ;
	}
	size = map.header.nx ;
	c = size / 2 ;
	fprintf(stderr, "Read model\n") ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rad = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		
		if (rad < rmin)
			map.f32_data[x*size*size + y*size + z] = -1.f ;
		else if (rad > rmax + 3)
			map.f32_data[x*size*size + y*size + z] = 0.f ;
		else if (rad > rmax - 3)
			map.f32_data[x*size*size + y*size + z] = -1.f ;
		else if (map.f32_data[x*size*size + y*size + z] < minmodel && map.f32_data[x*size*size + y*size + z] > -1000.)
			minmodel = map.f32_data[x*size*size + y*size + z] ;
	}
	fprintf(stderr, "Zeroed outer voxels and set inner voxels to be negative\n") ;
	fprintf(stderr, "Minimum value in annulus = %.3e\n", minmodel) ;
	
	sprintf(map.header.labels, "ResEx-phasing:zero_outer %s\n", extract_fname(argv[1])) ;
	map.header.rms = -1. ;
	sprintf(fname, "%s-zero.ccp4", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving to %s\n", fname) ;
	write_map(fname, &map) ;
	
	free_map(&map) ;
	
	return 0 ;
}
