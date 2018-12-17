#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "../../src/utils.h"
#include "../../src/volume.h"

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol ;
	float voxres, bfac, distsq ;
	float complex *model ;
	float *intens ;
	FILE *fp ;
	char fname[999] ;
	struct volume_data volume ;
	strcpy(volume.point_group, "222") ;
	
	if (argc < 4) {
		fprintf(stderr, "Sharpen: High pass filter by negative B-factor\n") ;
		fprintf(stderr, "----------------------------------------------\n") ;
		fprintf(stderr, "\nUsage: %s <cpx_fname> <voxres> <B-factor>\n", argv[0]) ;
		fprintf(stderr, "\twhere <voxres> is the resolution at 1 pixel in Angstroms\n") ;
		fprintf(stderr, "\tOptional: <point_group>. Currently support '222', '4', '1'\n") ;
		fprintf(stderr, "\nOutput: <cpx_fname>-sharp.cpx, <cpx_fname>-sharp-sym.raw\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float complex)) ;
	voxres = atof(argv[2]) ;
	bfac = atof(argv[3]) ;
	if (argc > 4)
		strcpy(volume.point_group, argv[4]) ;
	
	c = size / 2 ;
	vol = size*size*size ;
	// Adjust bfac to voxel^2 from A^2
	bfac /= (voxres*voxres) ;
	fprintf(stderr, "B-factor = %.3e vox^-2\n", bfac) ;
	// Adjust bfac for structure factors rather than intensities
	bfac /= 2. ;
	
	// Read model
	model = malloc(vol * sizeof(float complex)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// High pass filter
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		distsq = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c) ;
		// Divide by Debye-Waller factor
		model[x*size*size + y*size + z] *= expf(bfac * distsq) ;
	}
	
	// Write complex value to file
	sprintf(fname, "%s-sharp.cpx", remove_ext(argv[1])) ;
	fprintf(stderr, "Writing output to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Symmetrize intensities
	intens = malloc(vol * sizeof(float)) ;
	volume_symmetrize_centered(&volume, model, intens) ;
	
	// Write symmetrized intensities to file
	sprintf(fname, "%s-sharp-sym.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "Writing symmetrized intensity to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	free(intens) ;
	
	return 0 ;
}
