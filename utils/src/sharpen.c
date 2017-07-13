#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol ;
	float voxres, bfac, distsq ;
	float complex *model ;
	float *intens ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 5) {
		fprintf(stderr, "Format: %s <cpx_fname> <size> <voxres> <B-factor>\n", argv[0]) ;
		fprintf(stderr, "\twhere <voxres> is the resolution at 1 pixel in Angstroms\n") ;
		fprintf(stderr, "\tOptional: <point_group>. Currently support '222', '4', '1'\n") ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	voxres = atof(argv[3]) ;
	bfac = atof(argv[4]) ;
	
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
	if (argc > 5 && strcmp(argv[5], "4") == 0) {
		fprintf(stderr, "Symmetrizing by point group '4'\n") ;
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z)
			intens[x*size*size + y*size + z] = 0.25 * (powf(cabsf(model[x*size*size + y*size + z]), 2.f) +
													 powf(cabsf(model[x*size*size + (2*c-z)*size + y]), 2.f) +
													 powf(cabsf(model[x*size*size + (2*c-y)*size + (2*c-z)]), 2.f) +
													 powf(cabsf(model[x*size*size + z*size + (2*c-y)]), 2.f)) ;
	}
	else if (argc > 5 && strcmp(argv[5], "1") == 0) {
		fprintf(stderr, "Symmetrizing by point group '1'\n") ;
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z)
			intens[x*size*size + y*size + z] = powf(cabsf(model[x*size*size + y*size + z]), 2.f) ;
	}
	else {
		fprintf(stderr, "Symmetrizing by point group '222'\n") ;
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z)
			intens[x*size*size + y*size + z] = 0.25 * (powf(cabsf(model[x*size*size + y*size + z]), 2.f) +
													 powf(cabsf(model[x*size*size + y*size + (2*c-z)]), 2.f) +
													 powf(cabsf(model[x*size*size + (2*c-y)*size + z]), 2.f) +
													 powf(cabsf(model[x*size*size + (2*c-y)*size + (2*c-z)]), 2.f)) ;
	}
	
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
