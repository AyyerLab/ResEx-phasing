#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	long s, vol, mvol, x, y, z, vox, num_supp = 0 ;
	long mmin_x, mmax_x, ms_x ;
	long mmin_y, mmax_y, ms_y ;
	long mmin_z, mmax_z, ms_z ;
	float vox_x, vox_y, vox_z, *model, *mmodel ;
	float mean = 0.f, rms = 0.f, max = -1.e20, min = 1.e20 ;
	uint8_t *support ;
	struct ccp4_map map = {0} ;
	char fname[1024] ;
	FILE *fp ;
	
	// Parse command line arguments
	if (argc < 6) {
		fprintf(stderr, "Gen Map: Produce CCP4/MRC map from electron density\n") ;
		fprintf(stderr, "---------------------------------------------------\n") ;
		fprintf(stderr, "Needs voxel sizes of the density model\n") ;
		fprintf(stderr, "Also needs support mask file to calculat ebounding box\n") ;
		fprintf(stderr, "\nUsage: %s <model_fname> <vox_x> <vox_y> <vox_z> <supp_fname>\n", argv[0]) ;
		fprintf(stderr, "Set <supp_fname> = all if you want full volume in map\n") ;
		fprintf(stderr, "\nOutput: data/maps/<model_fname>.map.ccp4\n") ;
		return 1 ;
	}
	s = get_size(argv[1], sizeof(float)) ;
	vox_x = atof(argv[2]) ;
	vox_y = atof(argv[3]) ;
	vox_z = atof(argv[4]) ;
	vol = s*s*s ;
	
	// Read model and support
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	support = malloc(vol * sizeof(uint8_t)) ;
	if (strncmp(argv[5], "all", 3) == 0) {
		fprintf(stderr, "Assuming full cube to be within support\n") ;
		for (x = 0 ; x < vol ; ++x)
			support[x] = 1 ;
	}
	else {
		fp = fopen(argv[5], "rb") ;
		fread(support, sizeof(uint8_t), vol, fp) ;
		fclose(fp) ;
	}

	// From support, calculate size and position of bounding box
	// 	Also calculate model min, max and std inside support
	mmin_x = s ; mmin_y = s ; mmin_z = s ;
	mmax_x = 0 ; mmax_y = 0 ; mmax_z = 0 ;
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z) {
		vox = x*s*s + y*s + z ;
		if (support[vox] == 0)
			continue ;
		if (x < mmin_x)
			mmin_x = x ;
		if (y < mmin_y)
			mmin_y = y ;
		if (z < mmin_z)
			mmin_z = z ;
		if (x > mmax_x)
			mmax_x = x ;
		if (y > mmax_y)
			mmax_y = y ;
		if (z > mmax_z)
			mmax_z = z ;
		
		mean += model[vox] ;
		rms += model[vox] * model[vox] ;
		
		if (model[vox] > max)
			max = model[vox] ;
		if (model[vox] < min)
			min = model[vox] ;
		
		num_supp++ ;
	}
	free(support) ;
	
//	mmin_x -= 2 ; mmin_y -= 2 ; mmin_z -= 2 ;
//	mmax_x += 2 ; mmax_y += 2 ; mmax_z += 2 ;
	ms_x = mmax_x - mmin_x + 1 ;
	ms_y = mmax_y - mmin_y + 1 ;
	ms_z = mmax_z - mmin_z + 1 ;
	mvol = ms_x * ms_y * ms_z ;
	fprintf(stderr, "min = (%ld, %ld, %ld)\n", mmin_x, mmin_y, mmin_z) ;
	fprintf(stderr, "max = (%ld, %ld, %ld)\n", mmax_x, mmax_y, mmax_z) ;
	fprintf(stderr, "Model volume = %ld x %ld x %ld = %ld\n", ms_x, ms_y, ms_z, mvol) ;
	
	mean /= (float) num_supp ;
	rms /= (float) num_supp ;
	rms = sqrtf(rms - mean*mean) ;
	fprintf(stderr, "Effective solvent fraction = %f\n", 1. - ((double) num_supp) / mvol) ;
	fprintf(stderr, "max = %.6e, min = %.6e\n", max, min) ;
	fprintf(stderr, "mean = %.6e, rms = %.6e\n", mean, rms) ;
	
	// Extract sub-volume
	mmodel = malloc(mvol * sizeof(float)) ;
	for (x = 0 ; x < ms_x ; ++x)
	for (y = 0 ; y < ms_y ; ++y)
	for (z = 0 ; z < ms_z ; ++z)
//		mmodel[((z+ms_z/2)%ms_z)*ms_y*ms_x + ((y+ms_y/2)%ms_y)*ms_x + ((x+ms_x/2)%ms_x)] // ms_x*ms_y*ms_z, translated
		mmodel[z*ms_y*ms_x + y*ms_x + x] // ms_x*ms_y*ms_z
			= model[(x+mmin_x)*s*s + (y+mmin_y)*s + (z+mmin_z)] ; // s*s*s
	free(model) ;

	// Create map
	map.header.nx = ms_x ;
	map.header.ny = ms_y ;
	map.header.nz = ms_z ;
	map.header.mode = 2 ;
	map.header.nxstart = 0 ;
	map.header.nystart = 0 ;
	map.header.nzstart = 0 ;
	map.header.mx = ms_z ; // Note reversed axes
	map.header.my = ms_y ;
	map.header.mz = ms_x ;
	map.header.xlen = ms_z * vox_z ;
	map.header.ylen = ms_y * vox_y ;
	map.header.zlen = ms_x * vox_x ;
	map.header.alpha = 90.f ;
	map.header.beta = 90.f ;
	map.header.gamma = 90.f ;
	map.header.mapc = 3 ; // Due to reversed axes
	map.header.mapr = 2 ;
	map.header.maps = 1 ;
	map.header.dmax = max ;
	map.header.dmin = min ;
	map.header.dmean = mean ;
	map.header.ispg = 1 ;
	map.header.nsymbt = 0 ;
	memset(map.header.extra, 0, 100) ;
	map.header.xorig = 0.f ;
	map.header.yorig = 0.f ;
	map.header.zorig = 0.f ;
	strcpy(map.header.cmap, "MAP") ;
	map.header.cmap[3] = ' ' ;
	map.header.machst[0] = 0x44 ;
	map.header.machst[1] = 0x41 ; // Other two bytes zeroed
	map.header.rms = rms ;
	map.header.nlabl = 1 ;
	sprintf(map.header.labels[0], "ResEx-phasing:gen_map %s %ld %s", argv[1], s, argv[5]) ;
	map.data = malloc(mvol * sizeof(float)) ;
	memcpy(map.data, mmodel, mvol*sizeof(float)) ;
	
	// Write map to file
	fprintf(stderr, "box size = %f x %f x %f\n", map.header.xlen, map.header.ylen, map.header.zlen) ;
	sprintf(fname, "%s.ccp4", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving map file to %s\n", fname) ;
	write_map(fname, &map) ;
	
	// Free memory
	free(mmodel) ;
	free_map(&map) ;
	
	return 0 ;
}
