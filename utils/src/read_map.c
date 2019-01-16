#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	int nx, ny, nz, section_flag ;
	long psize, pvol, shx, shy, shz, x, y, z ;
	float px, py, pz, voxres[3] ;
	float central_sum = 0.f, edge_sum= 0.f ;
	float *padmodel ;
	struct ccp4_map map ;
	FILE *fp ;
	char fname[1024] ;
	
	// Parse command line arguments
	if (argc < 3) {
		fprintf(stderr, "Read Map: Parse CCP4/MRC map with given target voxel resolution\n") ;
		fprintf(stderr, "---------------------------------------------------------------\n") ;
		fprintf(stderr, "\nUsage: %s <map_fname> <voxres>\n", argv[0]) ;
		fprintf(stderr, "\twhere <voxres> is the resolution at 1 pixel in Angstroms\n") ;
		fprintf(stderr, "One can also give three different <voxres> parameters for different axes\n") ;
		fprintf(stderr, "\nOutput: data/<map_fname>-map.raw for just a dump of map\n") ;
		fprintf(stderr, "\tdata/convert/<map_fname>-<padded_size>.raw for padded model with approx. target voxel size\n") ;
		fprintf(stderr, "\tAlso output correction factors to get exact voxel size (to be used by 'fstretch')\n") ;
		return 1 ;
	}
	if (argc == 3) {
		voxres[0] = atoi(argv[2]) ;
		voxres[1] = atoi(argv[2]) ;
		voxres[2] = atoi(argv[2]) ;
	}
	else if (argc > 4) {
		voxres[0] = atoi(argv[2]) ;
		voxres[1] = atoi(argv[3]) ;
		voxres[2] = atoi(argv[4]) ;
	}
	else
		return 1 ;

	// Parse map file
	if (parse_map(argv[1], &map))
		return 1 ;

	if (map.header.mode != 2) {
		fprintf(stderr, "Need float32 map data\n") ;
		return 1 ;
	}
	nx = map.header.nx ;
	ny = map.header.ny ;
	nz = map.header.nz ;
	
	// Test if translation needed by comparing edge slice to central slice
	for (y = 0 ; y < ny ; ++y)
	for (z = 0 ; z < nz ; ++z) {
		edge_sum += fabsf(map.data[y*nz + z]) ;
		central_sum += fabsf(map.data[(nx/2)*ny*nz + y*nz + z]) ;
	}
	
	section_flag = edge_sum > central_sum ? 0 : 1 ;
	fprintf(stderr, "Edge vs central sums: (%e, %e)\n", edge_sum, central_sum) ;
	 
	// Calculate padded model size
	px = voxres[0]*map.header.mx/map.header.xlen ;
	py = voxres[1]*map.header.my/map.header.ylen ;
	pz = voxres[2]*map.header.mz/map.header.zlen ;
	px = px > py ? px : py ;
	px = px > pz ? px : pz ;
	psize = (int) px + 3 ;
	psize = psize%2 == 0 ? psize + 1 : psize ;
	pvol = psize*psize*psize ;
	
	fprintf(stderr, "Padded volume size = %ld\n", psize) ;
	px = voxres[0]*map.header.mx/map.header.xlen ;
	fprintf(stderr, "Ideal pad sizes = (%.3f, %.3f, %.3f)\n", pz, py, px) ;
	fprintf(stderr, "Stretch factors = (%.5f, %.5f, %.5f)\n", psize/pz, psize/py, psize/px) ;
	
	shx = (psize - nx) / 2 ;
	shy = (psize - ny) / 2 ;
	shz = (psize - nz) / 2 ;
	
	// Generate padded model
	padmodel = calloc(pvol, sizeof(float)) ;
	if (section_flag == 0) {
		fprintf(stderr, "Translating array by half its size\n") ;
		if (map.header.mapc == 1) {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = map.data[((x+nx/2)%nx)*ny*nz + ((y+ny/2)%ny)*nz + ((z+nz/2)%nz)] ;
		}
		else {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = map.data[((z+nz/2)%nz)*ny*nx + ((y+ny/2)%ny)*nx + ((x+nx/2)%nx)] ;
		}
	}
	else {
		fprintf(stderr, "Saving array without translation\n") ;
		if (map.header.mapc == 1) {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = map.data[x*ny*nz + y*nz + z] ;
		}
		else {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = map.data[z*ny*nx + y*nx + x] ;
		}
	}
	
	// Save to file
	sprintf(fname, "data/convert/%s-%ld.raw", remove_ext(extract_fname(argv[1])), psize) ;
	fprintf(stderr, "Saving padded model to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(padmodel, sizeof(float), pvol, fp) ;
	fclose(fp) ;
	
	free_map(&map) ;
	free(padmodel) ;
	
	return 0 ;
}
