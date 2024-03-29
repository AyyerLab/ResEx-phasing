#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long x, y, z, ix, iy, iz ;
	long size, fsize, cen, fcen, vol, fvol ;
	float tx, ty, tz, fx, fy, fz, cx, cy, cz, f[3] ;
	float complex *model, *fmodel ;
	float *temp ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 6) {
		fprintf(stderr, "Fourier Stretch: Slightly change q-sampling for given model\n") ;
		fprintf(stderr, "-----------------------------------------------------------\n") ;
		fprintf(stderr, "Uses linear interpolation assuming stretch factors <fx>, <fy>, <fz> close to 1\n") ;
		fprintf(stderr, "For larger rescaling use zero-padding in real-space\n") ;
		fprintf(stderr, "<fsize> is stretched model size\n") ;
		fprintf(stderr, "\nUsage: %s <cpx_fname> <fsize> <fx> <fy> <fz>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "\nOutput: <cpx_fname>_str_<fsize>.cpx (if <out_fname> is not given)\n") ;
		return 1 ;
	}
	
	size = get_size(argv[1], sizeof(float complex)) ;
	fsize = atoi(argv[2]) ;
	f[0] = atof(argv[3]) ;
	f[1] = atof(argv[4]) ;
	f[2] = atof(argv[5]) ;
	
	sprintf(fname, "%s_str_%ld.cpx", remove_ext(argv[1]), fsize) ;
	if (argc > 6)
		strcpy(fname, argv[6]) ;
	
	vol = size*size*size ;
	fvol = fsize*fsize*fsize ;
	cen = size / 2 ;
	fcen = fsize / 2 ;
	
	// Parse model
	model = malloc(vol * sizeof(float complex)) ;
	fmodel = malloc(fvol * sizeof(float complex)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed model\n") ;
	
	// For each voxel in fmodel, interpolate from model
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
			cx*cy*cz*model[ix*size*size + iy*size + iz] +
			cx*cy*fz*model[ix*size*size + iy*size + (iz+1)] +
			cx*fy*cz*model[ix*size*size + (iy+1)*size + iz] +
			cx*fy*fz*model[ix*size*size + (iy+1)*size + (iz+1)] +
			fx*cy*cz*model[(ix+1)*size*size + iy*size + iz] +
			fx*cy*fz*model[(ix+1)*size*size + iy*size + (iz+1)] +
			fx*fy*cz*model[(ix+1)*size*size + (iy+1)*size + iz] +
			fx*fy*fz*model[(ix+1)*size*size + (iy+1)*size + (iz+1)] ;
		
		if (y == fsize-1 && z == fsize-1)
			fprintf(stderr, "\rFinished x = %ld/%ld", x, fsize) ;
	}
	fprintf(stderr, "\n") ;
	
	// Write fmodel to file
	fprintf(stderr, "Writing to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(fmodel, sizeof(float complex), fvol, fp) ;
	fclose(fp) ;
	
	// Write symmetrized intensity to file
	temp = malloc(fvol * sizeof(float)) ;
	for (x = 0 ; x < fsize ; ++x)
	for (y = 0 ; y < fsize ; ++y)
	for (z = 0 ; z < fsize ; ++z)
		temp[x*fsize*fsize + y*fsize + z] = 0.25 * (pow(cabsf(fmodel[x*fsize*fsize + y*fsize + z]), 2.) + 
		                                         pow(cabsf(fmodel[x*fsize*fsize + y*fsize + (2*fcen-z)]), 2.) +
		                                         pow(cabsf(fmodel[x*fsize*fsize + (2*fcen-y)*fsize + z]), 2.) +
		                                         pow(cabsf(fmodel[x*fsize*fsize + (2*fcen-y)*fsize + (2*fcen-z)]), 2.)) ;
	
	sprintf(fname, "%s-sym.raw", remove_ext(fname)) ;
	fp = fopen(fname, "wb") ;
	fwrite(temp, sizeof(float), fvol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(model) ;
	free(fmodel) ;
	free(temp) ;
	
	return 0 ;
}
