#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char *argv[]) {
	long xs, ys, zs ;
	long i[3], x, y, z, num_points ;
	long center, size, spotnum[3], vol ;
	long t[3], a[3] ;
	double dist ;
	float *object ;
	char fname[100] ;
	FILE *fp ;
	float **cell_list ;
	int bin, num_bins, **count_list ;
	
	if (argc < 6) {
		fprintf(stderr, "Format: %s <model_fname> <size> <ax> <ay> <az>\n", argv[0]) ;
		return 1 ;
	}
	
	num_bins = 10 ;
	size = atoi(argv[2]) ; // Array size (=1501 for PS2)
	a[0] = atof(argv[3]) ; // Lattice constant in voxel units
	a[1] = atof(argv[4]) ; // Lattice constant in voxel units
	a[2] = atof(argv[5]) ; // Lattice constant in voxel units
	
	center = size / 2 ;
	vol = (long)size*size*size ;
	
	// Allocate cells
	xs = 2 * ((long) (a[0]/2)) + 1 ;
	ys = 2 * ((long) (a[1]/2)) + 1 ;
	zs = 2 * ((long) (a[2]/2)) + 1 ;
	cell_list = malloc(num_bins * sizeof(float*)) ;
	count_list = malloc(num_bins * sizeof(int*)) ;
	for (bin = 0 ; bin < num_bins ; ++bin) {
		cell_list[bin] = calloc(xs*ys*zs, sizeof(float)) ;
		count_list[bin] = calloc(xs*ys*zs, sizeof(int)) ;
	}
	
	fprintf(stderr, "Generated cells with size = %ldx%ldx%ld\n", xs, ys, zs) ;
	
	// Calculate number of spots to consider in each dimension
	for (x = 0 ; x < 3 ; ++x)
		spotnum[x] = (long) floor(center / a[x]) ;
	fprintf(stderr, "spotnums = %ld, %ld, %ld\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
	// Parse model
	object = calloc(vol, sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(object, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed model\n") ;
	
//	int sign ;
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
		
		// Determine peak location
		dist = 0. ;
		for (x = 0 ; x < 3 ; ++x) {
			t[x] = a[x] * i[x] ;
			dist += t[x] * t[x] ;
			t[x] += center ;
		}
		
		// Exclude non-Bragg region
		if (dist > 250. * 250.)
			continue ;
		// Exclude beamstop
		if (dist < 25. * 25.)
			continue ;
		
		// Determine bin
		bin = ((sqrt(dist) - 30.) * num_bins / (250. - 30.)) ;
		if (bin > num_bins - 1)
			continue ;
		
		num_points++ ;
		
		// Only for integer Bragg spacing
		for (x = 0 ; x < xs ; ++x)
		for (y = 0 ; y < ys ; ++y)
		for (z = 0 ; z < zs ; ++z) {
			cell_list[bin][x*ys*zs + y*zs + z] += 
				object[(int)(t[0]-xs/2+x)*size*size + (t[1]-ys/2+y)*size + (t[2]-zs/2+z)] ;
			count_list[bin][x*ys*zs + y*zs + z]++ ; 
		}
	}
	fprintf(stderr, "Generated cell using %ld peaks\n", num_points) ;
	
	// Write cells to files
	for (bin = 0 ; bin < num_bins ; ++bin) {
		for (x = 0 ; x < xs*ys*zs ; ++x)
		if (count_list[bin][x] > 0)
			cell_list[bin][x] /= count_list[bin][x] ;
		
		sprintf(fname, "data/cells/%ld_%ld_%ld_%.2d.raw", xs, ys, zs, bin) ;
		fp = fopen(fname, "wb") ;
		fwrite(cell_list[bin], sizeof(float), xs*ys*zs, fp) ;
		fclose(fp) ;
		
		free(cell_list[bin]) ;
		free(count_list[bin]) ;
	}
	
	free(object) ;
	free(cell_list) ;
	free(count_list) ;
	
	return 0 ;
}
