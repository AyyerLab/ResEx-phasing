#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

float gaussian(double x, float width) {
	return exp(- x*x / 2 / width/width) ;
}

void create_kernel(float *kernel, long ksize, float width) {
	long i, j, k, krad ;
	
	krad = (ksize - 1) / 2 ;
	
	for (i = 0 ; i < ksize ; ++i)
	for (j = 0 ; j < ksize ; ++j)
	for (k = 0 ; k < ksize ; ++k)
		kernel[i*ksize*ksize + j*ksize + k] = 
			gaussian(
				sqrt((krad-i)*(krad-i) + (krad-j)*(krad-j) + (krad-k)*(krad-k)), 
				width) ;
}

int main(int argc, char* argv[]) {
	long x, y, z, i, j, k ;
	long ksize, krad, size, vol, num_supp ;
	float blur, thresh ;
	float *density, *kernel ;
	uint8_t *support ;
	FILE *fp ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <density_fname> <blur> <threshold>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		return 1 ;
	}
	blur = atof(argv[2]) ;
	thresh = atof(argv[3]) ;
	
	size = 501 ;
	vol = size*size*size ;
	
	// Create Gaussian kernel of given width
	ksize = 3 * blur > 1 ? (long)3 * blur : 1 ;
	krad = (ksize - 1) / 2 ;
	kernel = malloc(ksize * ksize * ksize * sizeof(float)) ;
	create_kernel(kernel, ksize, blur) ;
	fprintf(stderr, "Generated kernel with ksize = %ld\n", ksize) ;
	
	// Parse density
	density = malloc(vol * sizeof(float)) ;
	support = calloc(vol, sizeof(uint8_t)) ;
	fp = fopen(argv[1], "rb") ;
	fread(density, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed density\n") ;
	
	// Convolve bright voxels with kernel, threshold and count voxels
	num_supp = 0 ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
	if (density[x*size*size + y*size + z] > 0.2 * thresh)
		for (i = 0 ; i < ksize ; ++i)
		for (j = 0 ; j < ksize ; ++j)
		for (k = 0 ; k < ksize ; ++k)
			if (density[x*size*size + y*size + z]*kernel[i*ksize*ksize + j*ksize + k] > thresh) {
				support[((i+x-krad)%size)*size*size + ((j+y-krad)%size)*size + ((k+z-krad)%size)] = 1 ;
				num_supp++ ;
			}
	fprintf(stderr, "Calculated support with %ld voxels\n", num_supp) ;
	
	// Write to file
	if (argc > 4)
		fp = fopen(argv[4], "wb") ;
	else
		fp = fopen("data/support_501.supp", "wb") ;
	fwrite(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(density) ;
	free(support) ;
	free(kernel) ;
	
	return 0 ;
}
