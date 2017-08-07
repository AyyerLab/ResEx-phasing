#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c ;
	float rad, rmin, rmax, minmodel = FLT_MAX ;
	float *model ;
	FILE *fp ;
	char fname[500] ;
	
	if (argc < 5) {
		fprintf(stderr, "Format: %s <model_fname> <size> <rmin> <rmax>\n", argv[0]) ;
		fprintf(stderr, "Optional: <subtract_min>. Default 0.\n") ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	rmin = atof(argv[3]) ;
	rmax = atof(argv[4]) ;
	
	c = size / 2 ;
	model = malloc(size*size*size * sizeof(float)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), size*size*size, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Read model\n") ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rad = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		
		if (rad < rmin)
			model[x*size*size + y*size + z] = -1.f ;
		else if (rad > rmax + 6)
			model[x*size*size + y*size + z] = 0.f ;
		else if (rad > rmax - 6)
			model[x*size*size + y*size + z] = -1.f ;
		else if (model[x*size*size + y*size + z] < minmodel && model[x*size*size + y*size + z] > -1000.)
			minmodel = model[x*size*size + y*size + z] ;
	}
	fprintf(stderr, "Zeroed outer voxels and set inner voxels to be negative\n") ;
	fprintf(stderr, "Minimum value in annulus = %.3e\n", minmodel) ;
	
	if (argc > 5 && atoi(argv[5]) > 0) {
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
	fwrite(model, sizeof(float), size*size*size, fp) ;
	fclose(fp) ;
	
	free(model) ;
	
	return 0 ;
}
