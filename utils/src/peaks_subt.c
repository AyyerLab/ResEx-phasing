#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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
	int bin, num_bins ;
	long i[3], x, y, z, num_points, vox ;
	long center, size, spotnum[3], vol ;
	long xs, ys, zs, t[3], a[3] ;
	double dist ;
	float ref, val, *object, **scale_list ;
	uint8_t *flag ;
	char fname[100] ;
	FILE *fp ;
	
	if (argc < 6) {
		fprintf(stderr, "Format: %s <model_fname> <size> <ax> <ay> <az>\n", argv[0]) ;
		return 1 ;
	}
	
	num_bins = 10 ;
	size = atoi(argv[2]) ; // Array size (= 501 for Anton's PS2)
	a[0] = atof(argv[3]) ; // Lattice constant in voxel units
	a[1] = atof(argv[4]) ; // Lattice constant in voxel units
	a[2] = atof(argv[5]) ; // Lattice constant in voxel units
	
	center = size / 2 ;
	vol = (long)size*size*size ;
	
	// Generate scale_list
	xs = 2 * ((long) (a[0]/2)) + 1 ;
	ys = 2 * ((long) (a[1]/2)) + 1 ;
	zs = 2 * ((long) (a[2]/2)) + 1 ;
	scale_list = malloc(num_bins * sizeof(float*)) ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		scale_list[bin] = malloc(xs*ys*zs * sizeof(float)) ;
	
	// Parse cells
	for (bin = 0 ; bin < num_bins ; ++bin) {
		sprintf(fname, "data/cells/%ld_%ld_%ld_%.2d.raw", xs, ys, zs, bin) ;
		fp = fopen(fname, "rb") ;
		fread(scale_list[bin], sizeof(float), xs*ys*zs, fp) ;
		fclose(fp) ;
		
/*		ref = scale_list[bin][0] ;
		
		// Subtract by corner voxel
		for (x = 0 ; x < xs*ys*zs ; ++x)
			scale_list[bin][x] -= ref ;
*/	}
	
	// Calculate number of spots to consider in each dimension
	for (x = 0 ; x < 3 ; ++x)
		spotnum[x] = (long) floor(center / a[x]) ;
	fprintf(stderr, "spotnums = %ld, %ld, %ld\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
	// Parse model
	object = malloc(vol * sizeof(float)) ;
	flag = calloc(vol, sizeof(uint8_t)) ;
	fp = fopen(argv[1], "rb") ;
	fread(object, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed model\n") ;
	
	// For each peak
	num_points = 0 ;
	for (i[0] = -spotnum[0] ; i[0] <= spotnum[0] ; ++i[0])
	for (i[1] = -spotnum[1] ; i[1] <= spotnum[1] ; ++i[1])
	for (i[2] = -spotnum[2] ; i[2] <= spotnum[2] ; ++i[2]) {
		// Skip if peak is forbidden
		if (i[0] == 0 && i[1] == 0 && i[2]%2 == 1)
			continue ;
		if (i[1] == 0 && i[2] == 0 && i[0]%2 == 1)
			continue ;
		if (i[2] == 0 && i[0] == 0 && i[1]%2 == 1)
			continue ;
		
		// Locate peak
		dist = 0. ;
		for (x = 0 ; x < 3 ; ++x) {
			t[x] = a[x] * i[x] ;
			dist += t[x] * t[x] ;
			t[x] += center ;
		}
		dist = sqrt(dist) ;
		
		// Exclude non-Bragg region
		if (dist > 250.)
			continue ;
		// Exclude beamstop
		if (dist < 25.)
			continue ;
		
		// Determine bin
		bin = ((dist - 30.) * num_bins / (250. - 30.)) ;
		if (bin > num_bins - 1)
			continue ;
		
		num_points++ ;
		
		// Subtract peak from intensity
		val = object[t[0]*size*size + t[1]*size + t[2]] 
		      / scale_list[bin][(xs/2)*ys*zs + (ys/2)*zs + (zs/2)] ;
		ref = scale_list[bin][0] ;
		
		for (x = 0 ; x < xs ; ++x)
		for (y = 0 ; y < ys ; ++y)
		for (z = 0 ; z < zs ; ++z) {
			vox = (int)(t[0]-xs/2+x)*size*size + (t[1]-ys/2+y)*size + (t[2]-zs/2+z) ;
			
			object[vox] -= val * (scale_list[bin][x*ys*zs + y*zs + z] - ref) ;
		}
		
		// Replace peak voxel by average of nearest neighbours
		vox = t[0]*size*size + t[1]*size + t[2] ;
		object[vox] = (object[vox+1] + object[vox-1]
		             + object[vox+size] + object[vox-size]
		             + object[vox+size*size] + object[vox-size*size]) / 6. ;
	}
	fprintf(stderr, "Subtracted %ld peaks\n", num_points) ;
	
	// Write damped model to file
	sprintf(fname, "%s_subt.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(object, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	for (bin = 0 ; bin < num_bins ; ++bin)
		free(scale_list[bin]) ;
	free(object) ;
	free(scale_list) ;
	free(flag) ;
	
	return 0 ;
}
