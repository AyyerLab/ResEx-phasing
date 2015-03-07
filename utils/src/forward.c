#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

float gaussian(float x, float width) {
	return exp(- x*x / 2 / width/width) ;
}

int main(int argc, char *argv[]) {
	long x, y, z, i[3], t[3], spotnum[3] ;
	long a[3], kxsize, kysize, kzsize ;
	long kxrad, kyrad, kzrad ;
	long size, c, vol, n_cell ;
	double dist ;
	float blur, val, b_factor ;
	float complex *model, cval ;
	float *intens, *kernel, *factor ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 6) {
		fprintf(stderr, "Format: %s <model_fname> <size> <ax> <ay> <az>\n", argv[0]) ;
		return 1 ;
	}
	
	size = atoi(argv[2]) ;
	a[0] = atoi(argv[3]) ;
	a[1] = atoi(argv[4]) ;
	a[2] = atoi(argv[5]) ;
	
	c = size / 2 ;
	vol = (long) size*size*size ;
	
	blur = 0.5 ; // Width of Bragg peak in voxels
	n_cell = 100 ; // Number of unit cells
//	b_factor = sqrt(log(n_cell) * 0.1 / M_PI / M_PI) ; // q-dependent factor
	b_factor = 0.8 ;
	fprintf(stderr, "n_cell = %ld, b_factor = %.4f\n", n_cell, b_factor) ;
	
	// Parse model of molecular transform
	model = malloc(vol * sizeof(float complex)) ;
	intens = calloc(vol, sizeof(float)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Create kernel
	kxsize = 2 * (a[0]/2) + 1 ;
	kysize = 2 * (a[1]/2) + 1 ;
	kzsize = 2 * (a[2]/2) + 1 ;
	kxrad = kxsize / 2 ;
	kyrad = kysize / 2 ;
	kzrad = kzsize / 2 ;
	
	kernel = malloc(kxsize * kysize * kzsize * sizeof(float)) ;
	float ktotal = 0.f ;
	
	for (x = 0 ; x < kxsize ; ++x)
	for (y = 0 ; y < kysize ; ++y)
	for (z = 0 ; z < kzsize ; ++z) {
		kernel[x*kysize*kzsize + y*kzsize + z] = 
			gaussian(sqrt((kxrad-x)*(kxrad-x) + 
			              (kyrad-y)*(kyrad-y) + 
			              (kzrad-z)*(kzrad-z)), 
			         blur) ;
		
		ktotal += kernel[x*kysize*kzsize + y*kzsize + z] ;
	}
	// Normalize kernel
	for (x = 0 ; x < kxsize*kysize*kzsize ; ++x)
		kernel[x] /= ktotal ;
	
	fprintf(stderr, "Created kernel with ksize: %ldx%ldx%ld\n", kxsize, kysize, kzsize) ;
	
	// Calculate factor
	factor = malloc(c * sizeof(float)) ;
	for (x = 0 ; x < c ; ++x)
//		factor[x] = exp(- 2. * M_PI * M_PI * b_factor * b_factor * x * x / c / c) ;
		factor[x] = exp(- M_PI * M_PI * b_factor * b_factor * x * x / c / c) ;
	
	// Calculate spotnums
	for (x = 0 ; x < 3 ; ++x)
		spotnum[x] = (long) floor(c / a[x]) ;
	fprintf(stderr, "spotnums = %ld, %ld, %ld\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
/*	// Geenerate list of hkl intensities with coherent sum
	// Convolve with Gaussian and add Wilson suppression
	for (i[0] = -spotnum[0] ; i[0] <= spotnum[0] ; ++i[0])
	for (i[1] = -spotnum[1] ; i[1] <= spotnum[1] ; ++i[1])
	for (i[2] = -spotnum[2] ; i[2] <= spotnum[2] ; ++i[2]) {
		// Skip if forbidden peak
		if (i[0] == 0 && i[1] == 0 && i[2]%2 == 1)
			continue ;
		if (i[1] == 0 && i[2] == 0 && i[0]%2 == 1)
			continue ;
		if (i[2] == 0 && i[0] == 0 && i[1]%2 == 1)
			continue ;
		
		// Determine location of pixel
		dist = 0. ;
		for (x = 0 ; x < 3 ; ++x) {
			t[x] = a[x] * i[x] ;
			dist += t[x] * t[x] ;
		}
		dist = sqrt(dist) ;
		
		if (dist > c - 6.)
			continue ;
		
		// Coherent sum
		cval = model[(c+t[0])*size*size + (c+t[1])*size + (c+t[2])] ;
		cval += model[(c-t[0])*size*size + (c-t[1])*size + (c+t[2])] * powf(-1.f, i[0]+i[1]) ;
		cval += model[(c-t[0])*size*size + (c+t[1])*size + (c-t[2])] * powf(-1.f, i[2]+i[0]) ;
		cval += model[(c+t[0])*size*size + (c-t[1])*size + (c-t[2])] * powf(-1.f, i[1]+i[2]) ;
		
		// Convert to intensity
		val = powf(cabsf(cval), 2.f) ;
		
		// Multiply by Wilson factor and number of unit cells
		val *= factor[(int)dist] * n_cell ;
		
		// Convolve with kernel and add to intens
		for (x = 0 ; x < kxsize ; ++x)
		for (y = 0 ; y < kysize ; ++y)
		for (z = 0 ; z < kzsize ; ++z)
			intens[(t[0]+c+x-kxrad) + (t[1]+c+y-kyrad)*size + (t[2]+c+z-kzrad)*size*size]
				+= val * kernel[x*kysize*kzsize + y*kzsize + z] ;
	}
*/	
	float eq, dq ;
	// Generate speckle intensities with incoherent sum
	// Apply Wilson factor and add to merge
	for (x = -c ; x < c ; ++x)
	for (y = -c ; y < c ; ++y)
	for (z = -c ; z < c ; ++z) {
		dist = sqrt(x*x + y*y + z*z) ;
		
		if (dist > c - 2.)
			continue ;
		
		// Incoherent sum
		val = powf(cabsf(model[(c+x)*size*size + (c+y)*size + (c+z)]), 2.f) ;
		val += powf(cabsf(model[(c-x)*size*size + (c-y)*size + (c+z)]), 2.f) ;
		val += powf(cabsf(model[(c-x)*size*size + (c+y)*size + (c-z)]), 2.f) ;
		val += powf(cabsf(model[(c+x)*size*size + (c-y)*size + (c-z)]), 2.f) ;
		
/*		// Multiply by Wilson factor
		val *= (1.f - factor[(int)dist]) ;
		
		// Add to intensity grid
		intens[(c+x) + (c+y)*size + (c+z)*size*size] += val ;
*/		
		eq = factor[(int) dist] ;
		
		dq = pow(1- eq*eq, 3.) ;
		dq /= 1 + eq*eq - 2*eq*cos(2*M_PI*x/a[0]) ;
		dq /= 1 + eq*eq - 2*eq*cos(2*M_PI*y/a[1]) ;
		dq /= 1 + eq*eq - 2*eq*cos(2*M_PI*z/a[2]) ;
		
		intens[(c+x) + (c+y)*size + (c+z)*size*size] += dq * val ;
	}
	
	// Write to file
	sprintf(fname, "%s_crap.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(model) ;
	free(intens) ;
	free(kernel) ;
	free(factor) ;
	
	return 0 ;
}
