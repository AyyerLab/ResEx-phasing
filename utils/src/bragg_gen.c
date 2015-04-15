#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>

int main(int argc, char *argv[]) {
	long i[3], t[3], spotnum[3], a[3] ;
	long x, size, c, vol, hkls[3], hklvol ;
	double dist ;
	float complex *model, *hkl ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 6) {
		fprintf(stderr, "Format: %s <model_fname> <size> <ax> <ay> <az>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		return 1 ;
	}
	
	// Needs (6,4,3) Bragg spacings
	size = atoi(argv[2]) ;
	a[0] = atoi(argv[3]) ;
	a[1] = atoi(argv[4]) ;
	a[2] = atoi(argv[5]) ;
	
	c = size / 2 ;
	vol = (long) size*size*size ;
	
	// Parse complex molecular transform
	model = malloc(vol * sizeof(float complex)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Calculate spotnums
	for (x = 0 ; x < 3 ; ++x)
		spotnum[x] = (long) floor(c / a[x]) ;
	fprintf(stderr, "spotnums = (%ld,%ld,%ld)\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
	// Allocate hkl array
	hklvol = 1 ;
	for (x = 0 ; x < 3 ; ++x) {
		hkls[x] = 2*spotnum[x] + 1 ;
		hklvol *= hkls[x] ;
	}
	fprintf(stderr, "hklsize =(%ld,%ld,%ld)\n", hkls[0], hkls[1], hkls[2]) ;
	
	hkl = malloc(hklvol * sizeof(float complex)) ;
	for (x = 0 ; x < hklvol ; ++x)
		hkl[x] = FLT_MAX ;
	
	// Sample peaks and save as reversed order
	for (i[0] = -spotnum[0] ; i[0] <= spotnum[0] ; ++i[0])
	for (i[1] = -spotnum[1] ; i[1] <= spotnum[1] ; ++i[1])
	for (i[2] = -spotnum[2] ; i[2] <= spotnum[2] ; ++i[2]) {
		// Determine location of pixel
		dist = 0. ;
		for (x = 0 ; x < 3 ; ++x) {
			t[x] = a[x] * i[x] ;
			dist += t[x] * t[x] ;
			t[x] += c ;
		}
		dist = sqrt(dist) ;
		
		if (dist > c - 1.)
			continue ;
		// Sample while reversing array
		hkl[(i[0]+spotnum[0]) + (i[1]+spotnum[1])*hkls[0] + (i[2]+spotnum[2])*hkls[1]*hkls[0]]
		 = model[t[0]*size*size + t[1]*size + t[2]] ;
	}
	
	// Save to file
	sprintf(fname, "data/hkl_%ld_%ld_%ld.cpx", hkls[2], hkls[1], hkls[0]) ;
	if (argc > 6)
		strcpy(fname, argv[6]) ;
	fp = fopen(fname, "wb") ;
	fwrite(hkl, sizeof(float complex), hklvol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(model) ;
	free(hkl) ;
	
	return 0 ;
}
