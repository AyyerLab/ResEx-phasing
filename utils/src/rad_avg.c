#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long x, y, z, vox, size, c, vol ;
	int num_bins, *count, bin ;
	double binsize, *rad_avg ;
	float *model ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <model_fname>\n", argv[0]) ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	
	vol = size*size*size ;
	c = size / 2 ;
	num_bins = ceil(sqrt(3.) * size / 2.) + 1 ;
	binsize = 1. ;
	
	model = malloc(vol * sizeof(float)) ;
	rad_avg = calloc(num_bins, sizeof(double)) ;
	count = calloc(num_bins, sizeof(int)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		vox = x*size*size + y*size + z ;
		bin = (int) (sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) / binsize) ;
		if (bin <  0 || bin >= num_bins) {
			fprintf(stderr, "bin[%ld, %ld, %ld] = %d > %d\n", x, y, z, bin, num_bins) ;
			continue ;
		}
		
		rad_avg[bin] += model[vox] ;
		count[bin]++ ;
	}
	
	for (x = 0 ; x < num_bins ; ++x)
	if (count[x] > 0)
		rad_avg[x] /= count[x] ;
	
//	for (x = 0 ; x < vol ; ++x)
//		model[x] = rad_avg[bin[x]] ;
	
	sprintf(fname, "%s-avg.bin", remove_ext(argv[1])) ;
	fprintf(stderr, "Writing radial average to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(rad_avg, sizeof(double), num_bins, fp) ;
	fclose(fp) ;
	
/*	sprintf(fname, "%s_avg.dat", remove_ext(argv[1])) ;
	fp = fopen(fname, "w") ;
	for (x = 0 ; x < num_bins ; ++x)
	if (count[x] > 0)
		fprintf(fp, "%ld\t%.6e\n", x, rad_avg[x]) ;
	fclose(fp) ;
*/
	free(model) ;
	free(rad_avg) ;
	free(count) ;
	
	return 0 ;
}
