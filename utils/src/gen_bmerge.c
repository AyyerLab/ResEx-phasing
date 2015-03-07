#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>

#define HMAX 39
#define KMAX 65
#define LMAX 88

int main(int argc, char *argv[]) {
	int bin, num_bins, vox ;
	int h, k, l, hsize, ksize, lsize, hklvol ;
	long x, y, z, num_points, center, size, vol ;
	long xs, ys, zs, t[3], a[3] ;
	double dist ;
	float val, *object, *intensarray, **scale_list ;
	double complex *hklarray ;
	uint8_t *flag ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 5) {
		fprintf(stderr, "Format: %s <hkl_fname> <ax> <ay> <az> <size>\n", argv[0]) ;
		return 1 ;
	}
	
	num_bins = 20 ;
	a[0] = atof(argv[2]) ; // Lattice constant in voxel units
	a[1] = atof(argv[3]) ; // Lattice constant in voxel units
	a[2] = atof(argv[4]) ; // Lattice constant in voxel units
	size = atoi(argv[5]) ; // Array size (= 501 for Anton's PS2)
	
	center = size / 2 ;
	vol = (long)size*size*size ;
	
	hsize = 2*HMAX + 1 ;
	ksize = 2*KMAX + 1 ;
	lsize = 2*LMAX + 1 ;
	hklvol = hsize * ksize * lsize ;
	
	// Allocate scale_list
	xs = 2 * ((long) (a[2]/2)) + 1 ;
	ys = 2 * ((long) (a[1]/2)) + 1 ;
	zs = 2 * ((long) (a[0]/2)) + 1 ;
	scale_list = malloc(num_bins * sizeof(float*)) ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		scale_list[bin] = malloc(xs*ys*zs * sizeof(float)) ;
	
	// Parse cells
	for (bin = 0 ; bin < num_bins ; ++bin) {
//		sprintf(fname, "../partial/data/cells/%ld_%ld_%ld_%.2d.raw", xs, ys, zs, bin) ;
		sprintf(fname, "../partial/data/cells/%ld_%ld_%ld_%.2d.raw", xs, ys, zs, 0) ;
		fp = fopen(fname, "rb") ;
		fread(scale_list[bin], sizeof(float), xs*ys*zs, fp) ;
		fclose(fp) ;
		
//		for (x = 0 ; x < xs*ys*zs ; ++x)
//			scale_list[bin][x] -= scale_list[bin][0] ;
	}
	
	// Parse complex hkl values
	hklarray = malloc(hklvol * sizeof(double complex)) ;
	fp = fopen(argv[1], "rb") ;
	fread(hklarray, sizeof(double complex), hklvol, fp) ;
	fclose(fp) ;
	
	// Calculate Bragg intensities and reverse array
	intensarray = malloc(hklvol * sizeof(float)) ;
	for (h = 0 ; h < hklvol ; ++h) 
		intensarray[h] = pow(cabs(hklarray[h]), 2.) ;
	
	fp = fopen("hkl_intens.raw", "wb") ;
	fwrite(intensarray, sizeof(float), hklvol, fp) ;
	fclose(fp) ;
	
	// Allocate 3D merge array
	object = calloc(vol, sizeof(float)) ;
	flag = calloc(vol, sizeof(uint8_t)) ;
	num_points = 0 ;
	
	// For each Bragg peak,
	for (h = 0 ; h < hsize ; ++h)
	for (k = 0 ; k < ksize ; ++k)
	for (l = 0 ; l < lsize ; ++l) {
		// Calculate lattice position
		t[0] = (h-HMAX) * a[0] ;
		t[1] = (k-KMAX) * a[1] ;
		t[2] = (l-LMAX) * a[2] ;
		
		dist = 0. ;
		for (x = 0 ; x < 3 ; ++x) {
			dist += t[x]*t[x] ;
			t[x] += center ;
		}
		dist = sqrt(dist) ;
		
		// Check for out-of-bounds
		if (t[0] < 0 || t[0] > size - 3)
			continue ;
		if (t[1] < 0 || t[1] > size - 3)
			continue ;
		if (t[2] < 0 || t[2] > size - 3)
			continue ;
		
		// Determine bin
		bin = ((dist - 30.) * num_bins / (250. - 30.)) ;
		if (bin > num_bins - 1 || bin < 0)
			continue ;
		
		num_points++ ;
		
		// Spread out Bragg peak using known shape
		// For overlaps, take mean of all contributions
		for (x = 0 ; x < xs ; ++x)
		for (y = 0 ; y < ys ; ++y)
		for (z = 0 ; z < zs ; ++z) {
/*		for (x = 1 ; x < xs-1 ; ++x)
		for (y = 1 ; y < ys-1 ; ++y)
		for (z = 1 ; z < zs-1 ; ++z) {
*/			vox = (int)(t[0]-xs/2+x)*size*size + (t[1]-ys/2+y)*size + (t[2]-zs/2+z) ;
			val = scale_list[bin][x*ys*zs + y*zs + z] * intensarray[h*ksize*lsize + k*lsize + l] ;
			
			if (flag[vox] == 0)
				object[vox] = val ;
			else
				object[vox] = (val + flag[vox]*object[vox]) / (flag[vox] + 1) ;
			
			flag[vox]++ ;
		}
	}
	
	fprintf(stderr, "Used %ld peaks to generate merge\n", num_points) ;
	
	// Reverse object indices
	float *tempobj = malloc(vol * sizeof(float)) ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		tempobj[x*size*size + y*size + z] = object[z*size*size + y*size + x] ;
	
	// Write object to file
	fp = fopen("data/bragg_merge.raw", "wb") ;
	fwrite(tempobj, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(object) ;
	free(tempobj) ;
	free(hklarray) ;
	free(intensarray) ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		free(scale_list[bin]) ;
	free(scale_list) ;
	
	return 0 ;
}
