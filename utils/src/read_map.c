#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	int nx, ny, nz, mode, mapx, mapy, mapz ;
	long psize, pvol, shx, shy, shz, x, y, z ;
	float lx, ly, lz, ax, ay, az ;
	float px, py, pz, voxres[3] ;
	float *model, *padmodel ;
	FILE *fp ;
	char fname[999] ;
	
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

	// Reading map
	// --------------------------------------------------------------------------------
	fp = fopen(argv[1], "rb") ;
	
	// Grid size
	fread(&nx, sizeof(int), 1, fp) ;
	fread(&ny, sizeof(int), 1, fp) ;
	fread(&nz, sizeof(int), 1, fp) ;
	fprintf(stderr, "Volume = (%d, %d, %d)\n", nx, ny, nz) ;
	
	// Mode
	fread(&mode, sizeof(int), 1, fp) ;
	if (mode == 0)
		fprintf(stderr, "Mode: (un)signed byte\n") ;
	else if (mode == 1)
		fprintf(stderr, "Mode: int16_t\n") ;
	else if (mode == 2)
		fprintf(stderr, "Mode: float32\n") ;
	else if (mode == 3)
		fprintf(stderr, "Mode: complex int16_t\n") ;
	else if (mode == 4)
		fprintf(stderr, "Mode: complex float32\n") ;
	else if (mode == 6)
		fprintf(stderr, "Mode: uint16_t\n") ;
	else if (mode == 16)
		fprintf(stderr, "Mode: RGB (3*uint8_t)\n") ;
	
//	fseek(fp, 12, SEEK_CUR) ;
	int mx, my, mz ;
	fread(&mx, sizeof(int), 1, fp) ;
	fread(&my, sizeof(int), 1, fp) ;
	fread(&mz, sizeof(int), 1, fp) ;
	fprintf(stderr, "Starting indices: (%d, %d, %d)\n", mx, my, mz) ;
	
	// Grid size
//	int mx, my, mz ;
	fread(&mx, sizeof(int), 1, fp) ;
	fread(&my, sizeof(int), 1, fp) ;
	fread(&mz, sizeof(int), 1, fp) ;
	fprintf(stderr, "Grid size = (%d, %d, %d)\n", mx, my, mz) ;
	
	// Cell size
	fread(&lx, sizeof(float), 1, fp) ;
	fread(&ly, sizeof(float), 1, fp) ;
	fread(&lz, sizeof(float), 1, fp) ;
	
	// Cell angles
	fread(&ax, sizeof(float), 1, fp) ;
	fread(&ay, sizeof(float), 1, fp) ;
	fread(&az, sizeof(float), 1, fp) ;
	fprintf(stderr, "Cell parameters: (%.3f, %.3f, %.3f) with angles (%.1f, %.1f, %.1f)\n",
	        lx, ly, lz, ax, ay, az) ;
	
	// Mapping of axes
	fread(&mapx, sizeof(int), 1, fp) ;
	fread(&mapy, sizeof(int), 1, fp) ;
	fread(&mapz, sizeof(int), 1, fp) ;
	fprintf(stderr, "xyz -> abc mapping: (%d, %d, %d)\n", mapx, mapy, mapz) ;
	
	// Value properties
	float min, max, mean ;
	fread(&min, sizeof(float), 1, fp) ;
	fread(&max, sizeof(float), 1, fp) ;
	fread(&mean, sizeof(float), 1, fp) ;
	
	fseek(fp, 128, SEEK_CUR) ;
	float rms ;
	fread(&rms, sizeof(float), 1, fp) ;
	fprintf(stderr, "min, max, mean, rms = (%.3f, %.3f, %.3f, %.3f)\n", min, max, mean, rms) ;
	fseek(fp, 804, SEEK_CUR) ;
	
	// Voxels
	model = malloc(nx*ny*nz* sizeof(float)) ;
	fread(model, sizeof(float), nx*ny*nz, fp) ;
	fclose(fp) ;
	// --------------------------------------------------------------------------------

	// Comparing edge slices to central slices
	float central_sum = 0.f, edge_sum= 0.f ;
	
	for (y = 0 ; y < ny ; ++y)
	for (z = 0 ; z < nz ; ++z) {
		edge_sum += fabsf(model[y*nz + z]) ;
		central_sum += fabsf(model[(nx/2)*ny*nz + y*nz + z]) ;
	}
	
	int section_flag = edge_sum > central_sum ? 0 : 1 ;
	fprintf(stderr, "sums: (%e, %e)\n", edge_sum, central_sum) ;
	// --------------------------------------------------------------------------------
	 
	/*
	sprintf(fname, "data/%s-map.raw", remove_ext(extract_fname(argv[1]))) ;
	fprintf(stderr, "Saving model to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), nx*ny*nz, fp) ;
	fclose(fp) ;
	*/
	
	px = voxres[0]*mx/lx ;
	py = voxres[1]*my/ly ;
	pz = voxres[2]*mz/lz ;
	px = px > py ? px : py ;
	px = px > pz ? px : pz ;
	psize = (int) px + 3 ;
	psize = psize%2 == 0 ? psize + 1 : psize ;
	pvol = psize*psize*psize ;
	
	fprintf(stderr, "Padded volume size = %ld\n", psize) ;
	px = voxres[0]*mx/lx ;
	fprintf(stderr, "Ideal pad sizes = (%.3f, %.3f, %.3f)\n", pz, py, px) ;
	fprintf(stderr, "Stretch factors = (%.5f, %.5f, %.5f)\n", psize/pz, psize/py, psize/px) ;
	
	shx = (psize - nx) / 2 ;
	shy = (psize - ny) / 2 ;
	shz = (psize - nz) / 2 ;
	
	padmodel = calloc(pvol, sizeof(float)) ;
	if (section_flag == 0) {
		fprintf(stderr, "Rotating array by half its size\n") ;
		if (mapx == 1) {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = model[((x+nx/2)%nx)*ny*nz + ((y+ny/2)%ny)*nz + ((z+nz/2)%nz)] ;
		}
		else {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = model[((z+nz/2)%nz)*ny*nx + ((y+ny/2)%ny)*nx + ((x+nx/2)%nx)] ;
		}
	}
	else {
		fprintf(stderr, "Saving array without rotation\n") ;
		if (mapx == 1) {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = model[x*ny*nz + y*nz + z] ;
		}
		else {
			for (x = 0 ; x < nx ; ++x)
			for (y = 0 ; y < ny ; ++y)
			for (z = 0 ; z < nz ; ++z)
				padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
				 = model[z*ny*nx + y*nx + x] ;
		}
	}
	
	sprintf(fname, "data/convert/%s-%ld.raw", remove_ext(extract_fname(argv[1])), psize) ;
	fprintf(stderr, "Saving padded model to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(padmodel, sizeof(float), pvol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	free(padmodel) ;
	
	return 0 ;
}
