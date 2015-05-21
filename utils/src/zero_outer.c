#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
	long x, y, z, size, c, rsq, rmax, rmaxsq ;
	float *model ;
	FILE *fp ;
	char fname[500] ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <model_fname> <size> <rmax>\n", argv[0]) ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	rmax = atoi(argv[3]) ;
	
	c = size / 2 ;
	rmaxsq = rmax*rmax ;
	
	model = malloc(size*size*size * sizeof(float)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), size*size*size, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Read model\n") ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rsq = (x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c) ;
		
		if (rsq < rmaxsq)
			model[x*size*size + y*size + z] = -1.f ;
		if (rsq > 0.95*c*c)
			model[x*size*size + y*size + z] = 0.f ;
	}
	fprintf(stderr, "Zeroed outer voxels and set inner voxels to be negative\n") ;
	
	sprintf(fname, "%s_zero.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "%s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), size*size*size, fp) ;
	fclose(fp) ;
	
	free(model) ;
	
	return 0 ;
}
		
