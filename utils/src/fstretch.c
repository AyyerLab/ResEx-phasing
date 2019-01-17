#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	long x, y, z, ix, iy, iz ;
	long size, fsize, cen, fcen, fvol ;
	float tx, ty, tz, fx, fy, fz, cx, cy, cz, f[3] ;
	float complex *fmodel ;
	float *temp ;
	char fname[1024], label[800] ;
	struct ccp4_map map = {0} ;
	
	if (argc < 6) {
		fprintf(stderr, "Fourier Stretch: Slightly change q-sampling for given model\n") ;
		fprintf(stderr, "-----------------------------------------------------------\n") ;
		fprintf(stderr, "Uses linear interpolation assuming stretch factors <fx>, <fy>, <fz> close to 1\n") ;
		fprintf(stderr, "For larger rescaling use zero-padding in real-space\n") ;
		fprintf(stderr, "<fsize> is stretched model size\n") ;
		fprintf(stderr, "\nUsage: %s <fmodel_name> <fsize> <fx> <fy> <fz>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "\nOutput: <fmodel_name>_str_<fsize>.ccp4 (if <out_fname> is not given)\n") ;
		return 1 ;
	}
	
	fsize = atoi(argv[2]) ;
	f[0] = atof(argv[3]) ;
	f[1] = atof(argv[4]) ;
	f[2] = atof(argv[5]) ;
	
	sprintf(fname, "%s_str_%ld.ccp4", remove_ext(argv[1]), fsize) ;
	if (argc > 6)
		strcpy(fname, argv[6]) ;
	
	fvol = fsize*fsize*fsize ;
	fcen = fsize / 2 ;
	
	// Parse model
	fmodel = malloc(fvol * sizeof(float complex)) ;
	
	parse_map(argv[1], &map) ;
	if (map.c64_data == NULL) {
		fprintf(stderr, "Need float complex data in input fmodel\n") ;
		return 1 ;
	}
	if (map.header.nx != map.header.ny || map.header.nx != map.header.nz) {
		fprintf(stderr, "Need cubic volume of data (input shape: %d, %d, %d)\n", map.header.nx, map.header.ny, map.header.nz) ;
		return 1 ;
	}
	size = map.header.nx ;
	cen = size / 2 ;
	fprintf(stderr, "Parsed model\n") ;
	
	// For each voxel in fmodel, interpolate from map
	// Keep center common
	for (x = 0 ; x < fsize ; ++x)
	for (y = 0 ; y < fsize ; ++y)
	for (z = 0 ; z < fsize ; ++z) {
		// Find position of voxel in model
		tx = (x-fcen)*f[0] + cen ;
		ty = (y-fcen)*f[1] + cen ;
		tz = (z-fcen)*f[2] + cen ;
		
		// Check if out of bounds
		if (tx < 0 || ty < 0 || tz < 0)
			continue ;
		if (tx > size-2 || ty > size-2 || tz > size-2)
			continue ;
		
		// Calculate interpolation weights
		ix = tx ;
		fx = tx - ix ;
		cx = 1. - fx ;
		iy = ty ;
		fy = ty - iy ;
		cy = 1. - fy ;
		iz = tz ;
		fz = tz - iz ;
		cz = 1. - fz ;
		
//		fmodel[z*fsize*fsize + y*fsize + x] = // To reverse indices
		fmodel[x*fsize*fsize + y*fsize + z] = 
			cx*cy*cz*map.c64_data[ix*size*size + iy*size + iz] +
			cx*cy*fz*map.c64_data[ix*size*size + iy*size + (iz+1)] +
			cx*fy*cz*map.c64_data[ix*size*size + (iy+1)*size + iz] +
			cx*fy*fz*map.c64_data[ix*size*size + (iy+1)*size + (iz+1)] +
			fx*cy*cz*map.c64_data[(ix+1)*size*size + iy*size + iz] +
			fx*cy*fz*map.c64_data[(ix+1)*size*size + iy*size + (iz+1)] +
			fx*fy*cz*map.c64_data[(ix+1)*size*size + (iy+1)*size + iz] +
			fx*fy*fz*map.c64_data[(ix+1)*size*size + (iy+1)*size + (iz+1)] ;
		
		if (y == fsize-1 && z == fsize-1)
			fprintf(stderr, "\rFinished x = %ld/%ld", x, fsize) ;
	}
	fprintf(stderr, "\n") ;
	
	// Write fmodel to file
	fprintf(stderr, "Writing to %s\n", fname) ;
	int sizes[3] = {fsize, fsize, fsize} ;
	float vsizes[3] = {-1.f, -1.f, -1.f} ; // Negative to indicate Fourier space
	sprintf(label, "ResEx-phasing:fstretch %s %ld %f %f %f\n", extract_fname(argv[1]), fsize, f[0], f[1], f[2]) ;
	int flipped = (map.header.mapc == 3) ;
	save_cpx_as_map(fname, fmodel, sizes, vsizes, label, flipped) ;
	
	// Write symmetrized intensity to file
	temp = malloc(fvol * sizeof(float)) ;
	for (x = 0 ; x < fsize ; ++x)
	for (y = 0 ; y < fsize ; ++y)
	for (z = 0 ; z < fsize ; ++z)
		temp[x*fsize*fsize + y*fsize + z] = 0.25 * (pow(cabsf(fmodel[x*fsize*fsize + y*fsize + z]), 2.) + 
		                                         pow(cabsf(fmodel[x*fsize*fsize + y*fsize + (2*fcen-z)]), 2.) +
		                                         pow(cabsf(fmodel[x*fsize*fsize + (2*fcen-y)*fsize + z]), 2.) +
		                                         pow(cabsf(fmodel[x*fsize*fsize + (2*fcen-y)*fsize + (2*fcen-z)]), 2.)) ;
	
	sprintf(fname, "%s-sym.ccp4", remove_ext(fname)) ;
	sprintf(label, "ResEx-phasing:fstretch(sym) %s %ld %f %f %f\n", extract_fname(argv[1]), fsize, f[0], f[1], f[2]) ;
	save_vol_as_map(fname, temp, sizes, vsizes, label, flipped) ;
	
	// Free memory
	free_map(&map) ;
	free(fmodel) ;
	free(temp) ;
	
	return 0 ;
}
