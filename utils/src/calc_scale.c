#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, num_vox = 0 ;
	double dot = 0., normsq = 0., rsq, scale ;
	float *obs_mag, *model_mag ;
	float complex *model ;
	FILE *fp ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <merge_fname>\n", argv[0]) ;
		return 1 ;
	}
	
	size = 501 ;
	c = size / 2 ;
	vol = size*size*size ;
	
	// Parse complex model
	model = malloc(vol * sizeof(float complex)) ;
	fp = fopen("data/4pbu_str_501.cpx", "rb") ;
	fread(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Calculate magnitude
	model_mag = malloc(vol * sizeof(float)) ;
	for (x = 0 ; x < vol ; ++x)
		model_mag[x] = cabsf(model[x]) ;
	
	free(model) ;
	
	// Parse experimental merge
	obs_mag = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(obs_mag, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate magnitude
	for (x = 0 ; x < vol ; ++x)
	if (obs_mag[x] > 0.)
		obs_mag[x] = sqrt(obs_mag[x]) ;
	
	// Calculate scale factor
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rsq = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c) ;
		
		if (rsq > 120.*120. && rsq < 150.*150. && obs_mag[x*size*size + y*size + z] > 0.) {
			dot += obs_mag[x*size*size + y*size + z] * model_mag[x*size*size + y*size + z] ;
			normsq += pow(obs_mag[x*size*size + y*size + z], 2.) ;
			num_vox++ ;
		}
	}
	
	scale = dot / normsq ;
	
	printf("Scale factor for magnitude: %.4e\n", scale) ;
	printf("%'ld voxels contributed to this merge\n", num_vox) ;
	
	// Free memory
	free(model_mag) ;
	free(obs_mag) ;
	
	return 0 ;
}
