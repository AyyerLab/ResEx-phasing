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
	int8_t *support ;
	char fname[1024], label[800] ;
	FILE *fp ;
	
	// Parse command line arguments
	if (argc < 6) {
		fprintf(stderr, "Gen Map: Produce CCP4/MRC map from electron density\n") ;
		fprintf(stderr, "---------------------------------------------------\n") ;
		fprintf(stderr, "Needs voxel sizes of the density model\n") ;
		fprintf(stderr, "Also needs support mask file to calculat ebounding box\n") ;
		fprintf(stderr, "\nUsage: %s <model_fname> <vox_x> <vox_y> <vox_z> <supp_fname>\n", argv[0]) ;
		fprintf(stderr, "Set <supp_fname> = all if you want full volume in map\n") ;
		fprintf(stderr, "\nOutput: <model_fname>.ccp4\n") ;
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
	
	support = malloc(vol * sizeof(int8_t)) ;
	if (strncmp(argv[5], "all", 3) == 0) {
		fprintf(stderr, "Assuming full cube to be within support\n") ;
		for (x = 0 ; x < vol ; ++x)
			support[x] = 1 ;
	}
	else {
		fp = fopen(argv[5], "rb") ;
		fread(support, sizeof(int8_t), vol, fp) ;
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
		
		num_supp++ ;
	}
	
//	mmin_x -= 2 ; mmin_y -= 2 ; mmin_z -= 2 ;
//	mmax_x += 2 ; mmax_y += 2 ; mmax_z += 2 ;
	ms_x = mmax_x - mmin_x + 1 ;
	ms_y = mmax_y - mmin_y + 1 ;
	ms_z = mmax_z - mmin_z + 1 ;
	mvol = ms_x * ms_y * ms_z ;
	fprintf(stderr, "min = (%ld, %ld, %ld)\n", mmin_x, mmin_y, mmin_z) ;
	fprintf(stderr, "max = (%ld, %ld, %ld)\n", mmax_x, mmax_y, mmax_z) ;
	fprintf(stderr, "Model volume = %ld x %ld x %ld = %ld\n", ms_x, ms_y, ms_z, mvol) ;
	
	fprintf(stderr, "Effective solvent fraction = %f\n", 1. - ((double) num_supp) / mvol) ;
	
	// Extract sub-volume
	mmodel = malloc(mvol * sizeof(float)) ;
	for (x = 0 ; x < ms_x ; ++x)
	for (y = 0 ; y < ms_y ; ++y)
	for (z = 0 ; z < ms_z ; ++z)
//		mmodel[((z+ms_z/2)%ms_z)*ms_y*ms_x + ((y+ms_y/2)%ms_y)*ms_x + ((x+ms_x/2)%ms_x)] // ms_x*ms_y*ms_z, translated
		mmodel[z*ms_y*ms_x + y*ms_x + x] // ms_x*ms_y*ms_z
			= model[(x+mmin_x)*s*s + (y+mmin_y)*s + (z+mmin_z)] ; // s*s*s
	
	// Write map to file
	sprintf(label, "ResEx-phasing:gen_map %s %ld %s", argv[1], s, argv[5]) ;
	sprintf(fname, "%s.ccp4", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving map file to %s\n", fname) ;
	int size[3] = {ms_x, ms_y, ms_z} ;
	float vsize[3] = {vox_x, vox_y, vox_z} ;
	save_vol_as_map(fname, mmodel, size, vsize, label, 1) ;
	
	// Free memory
	free(support) ;
	free(model) ;
	free(mmodel) ;
	
	return 0 ;
}
