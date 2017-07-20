#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <math.h>
#include <locale.h>

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, bin, num_bins, vox, *num_vox ;
	double *dot, *normsq, rad, binsize ;
	float *obs_mag, *model_mag ;
	FILE *fp ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <sym_model_fname> <merge_fname> <size>\n", argv[0]) ;
		return 1 ;
	}
	size = atoi(argv[3]) ;
	binsize = 2. ;
	
	setlocale(LC_ALL, "C") ; // For commas in large integers
	c = size / 2 ;
	vol = size*size*size ;
	
	// Parse complex model
	model_mag = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model_mag, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate magnitude
	for (x = 0 ; x < vol ; ++x)
		model_mag[x] = sqrt(model_mag[x]) ;
	
	// Parse experimental merge
	obs_mag = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[2], "rb") ;
	fread(obs_mag, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate magnitude
	for (x = 0 ; x < vol ; ++x)
	if (obs_mag[x] > 0.)
		obs_mag[x] = sqrt(obs_mag[x]) ;
	
	// Calculate scale factor
	double cx, cy, cz ;
	num_bins = (int) (size / binsize) ;
	dot = calloc(3*num_bins, sizeof(double)) ;
	normsq = calloc(3*num_bins, sizeof(double)) ;
	num_vox = calloc(3*num_bins, sizeof(long)) ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rad = sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		bin = (int) (rad / binsize) ;
		cx = fabs((x-c)) / rad ;
		cy = fabs((y-c)) / rad ;
		cz = fabs((z-c)) / rad ;
		vox = x*size*size + y*size + z ;
		if (obs_mag[vox] > 0. && cx > 0.95) {
			dot[0*num_bins + bin] += obs_mag[vox] * model_mag[vox] ;
			normsq[0*num_bins + bin] += pow(obs_mag[vox], 2.) ;
			num_vox[0*num_bins + bin]++ ;
		}
		else if (obs_mag[vox] > 0. && cy > 0.95) {
			dot[1*num_bins + bin] += obs_mag[vox] * model_mag[vox] ;
			normsq[1*num_bins + bin] += pow(obs_mag[vox], 2.) ;
			num_vox[1*num_bins + bin]++ ;
		}
		else if (obs_mag[vox] > 0. && cz > 0.95) {
			dot[2*num_bins + bin] += obs_mag[vox] * model_mag[vox] ;
			normsq[2*num_bins + bin] += pow(obs_mag[vox], 2.) ;
			num_vox[2*num_bins + bin]++ ;
		}
		if (y == size-1 && z == size-1)
			fprintf(stderr, "\rFinished x = %ld", x+1) ;
	}
	fprintf(stderr, "\n") ;
	
	// Write output to file
	fp = fopen("data/scale_q.dat", "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin) {
		if (normsq[0*num_bins + bin] > 0.)
			dot[0*num_bins + bin] /= normsq[0*num_bins + bin] ;
		if (normsq[1*num_bins + bin] > 1.)
			dot[1*num_bins + bin] /= normsq[1*num_bins + bin] ;
		if (normsq[2*num_bins + bin] > 2.)
			dot[2*num_bins + bin] /= normsq[2*num_bins + bin] ;
		fprintf(fp, "%8.3f %.6e %8ld %.6e %8ld %.6e %8ld\n", bin*binsize, 
				dot[0*num_bins + bin], num_vox[0*num_bins + bin],
				dot[1*num_bins + bin], num_vox[1*num_bins + bin],
				dot[2*num_bins + bin], num_vox[2*num_bins + bin]) ;
	}
	
	// Free memory
	free(model_mag) ;
	free(obs_mag) ;
	free(dot) ;
	free(normsq) ;
	free(num_vox) ;
	
	return 0 ;
}
