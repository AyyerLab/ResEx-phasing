#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <complex.h>

int main(int argc, char *argv[]) {
	long i[3], it[2][3], spotnum[3] ;
	long x, y, z, size, c, vol, hkls[3], hklvol ;
	float t[3], a[3], f[2][3], dist ;
	float complex *model, *hkl, val ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 6) {
		fprintf(stderr, "Format: %s <model_fname> <size> <ax> <ay> <az>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		return 1 ;
	}
	
	// For dominik-2.raw, (2.229, 3.055, 5.308)
	size = atoi(argv[2]) ;
	a[0] = atof(argv[3]) ;
	a[1] = atof(argv[4]) ;
	a[2] = atof(argv[5]) ;
	fprintf(stderr, "Parsed a = (%.3f, %.3f, %.3f)\n", a[0], a[1], a[2]) ;
	
	c = size / 2 ;
	vol = (long) size*size*size ;
	
	// Parse complex molecular transform
	model = malloc(vol * sizeof(float complex)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Calculate spotnums
	// For int_3x merge, spotnums = (83, 62, 41)
	spotnum[0] = 83 ;
	spotnum[1] = 62 ;
	spotnum[2] = 41 ;
//	for (x = 0 ; x < 3 ; ++x)
//		spotnum[x] = (long) floor(c / a[x]) ;
	fprintf(stderr, "spotnums = (%ld, %ld, %ld)\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
	// Allocate hkl array
	hklvol = 1 ;
	for (x = 0 ; x < 3 ; ++x) {
		hkls[x] = 2*spotnum[x] + 1 ;
		hklvol *= hkls[x] ;
	}
	fprintf(stderr, "hklsize =(%ld, %ld, %ld)\n", hkls[0], hkls[1], hkls[2]) ;
	
	hkl = malloc(hklvol * sizeof(float complex)) ;
	for (x = 0 ; x < hklvol ; ++x)
		hkl[x] = FLT_MAX ;
	
	// Sample peaks and save as reversed order
	for (i[0] = 0 ; i[0] < hkls[0] ; ++i[0])
	for (i[1] = 0 ; i[1] < hkls[1] ; ++i[1])
	for (i[2] = 0 ; i[2] < hkls[2] ; ++i[2]) {
		// Determine location of pixel
		dist = 0. ;
		for (x = 0 ; x < 3 ; ++x) {
			t[x] = a[x] * (i[x] - spotnum[x]) ;
			dist += t[x] * t[x] ;
			t[x] += c ;
		}
		dist = sqrtf(dist) ;
		
//		fprintf(stderr, "t = (%.3f, %.3f, %.3f) ", t[0], t[1], t[2]) ;
		
		if (dist > 220.)
			continue ;
		
		// Calculate interpolation factors
		for (x = 0 ; x < 3 ; ++x) {
			it[0][x] = (long) t[x] ;
			it[1][x] = it[0][x] + 1 ;
			f[1][x] = t[x] - it[0][x] ;
			f[0][x] = 1. - f[1][x] ;
		}
		
//		fprintf(stderr, "it[0] = (%ld, %ld, %ld)\n", it[0][0], it[0][1], it[0][2]) ;
//		fprintf(stderr, "f[0] = (%.3f, %.3f, %.3f)\n", f[0][0], f[0][1], f[0][2]) ;
		
		// Interpolate about t[x]
		val = 0.f ;
		for (x = 0 ; x < 2 ; ++x) 
		for (y = 0 ; y < 2 ; ++y) 
		for (z = 0 ; z < 2 ; ++z)
			val += f[x][0] * f[y][1] * f[z][2] * model[it[x][0]*size*size + it[y][1]*size + it[z][2]] ;
		
		hkl[i[0]*hkls[1]*hkls[2] + i[1]*hkls[2] + i[2]] = val ;
	}
	
	// Save to file
	sprintf(fname, "data/hkl_%ld_%ld_%ld.cpx", hkls[0], hkls[1], hkls[2]) ;
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
