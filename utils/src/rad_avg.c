#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
	int num_bins, *count, *bin ;
	float binsize, *model, *rad_avg ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <model_fname> <size>\n", argv[0]) ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	
	vol = size*size*size ;
	c = size / 2 ;
	num_bins = size ;
	binsize = 1. ;
	
	model = malloc(vol * sizeof(float)) ;
	bin = malloc(vol * sizeof(int)) ;
	rad_avg = calloc(num_bins, sizeof(float)) ;
	count = calloc(num_bins, sizeof(int)) ;
	
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		bin[x*size*size + y*size + z]
		 = (int) (sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) / binsize) ;
	
	for (x = 0 ; x < vol ; ++x) {
		if (bin[x] >= num_bins)
			continue ;
		
		rad_avg[bin[x]] += model[x] ;
		count[bin[x]]++ ;
	}
	
	for (x = 0 ; x < num_bins ; ++x)
	if (count[x] > 0)
		rad_avg[x] /= count[x] ;
	
	for (x = 0 ; x < vol ; ++x)
		model[x] = rad_avg[bin[x]] ;
	
/*	sprintf(fname, "%s_avg.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
*/	
	sprintf(fname, "%s_avg.dat", remove_ext(argv[1])) ;
	fp = fopen(fname, "w") ;
	for (x = 0 ; x < num_bins ; ++x)
	if (count[x] > 0)
		fprintf(fp, "%ld\t%.6e\n", x, rad_avg[x]) ;
	fclose(fp) ;
	
	free(model) ;
	free(rad_avg) ;
	free(count) ;
	free(bin) ;
	
	return 0 ;
}
