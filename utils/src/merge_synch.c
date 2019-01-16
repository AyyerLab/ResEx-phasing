#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <stdint.h>
#include "../../src/utils.h"

void make_rot(double phi, double rot[3][3]) {
	double c = cos(phi) ;
	double s = sin(phi) ;
	
	rot[0][0] = c ;
	rot[0][1] = 0 ;
	rot[0][2] = s ;
	rot[1][0] = 0 ;
	rot[1][1] = 1 ;
	rot[1][2] = 0 ;
	rot[2][0] = -s ;
	rot[2][1] = 0 ;
	rot[2][2] = c ;
}

int main(int argc, char *argv[]) {
	int i, j, det_flag, num_files = 0, avecount ;
	long x, y, t, fsize, size, center, vox  ;
	double cenx, ceny, pixsize, detd, lambda, phi, mscale, factor ;
	double rot[3][3], *det = NULL ;
	float *model3d, *weight, ave ;
	unsigned short *frame ;
	char fname[1024], header[512] ;
	FILE *fp, *fp_list ;
	long z ;
	double tx, ty, tz, fx, fy, fz, cx, cy, cz, w, f ;
	double rot_pix[3] ;
	int8_t *mask ;
	
	if (argc < 2) {
		fprintf(stderr, "Merge Synch: Basic frame merging using CBF file list\n") ;
		fprintf(stderr, "----------------------------------------------------\n") ;
		fprintf(stderr, "Applies radial average subtraction and merges in lab frame\n") ;
		fprintf(stderr, "File list is text file with one path per line\n") ;
		fprintf(stderr, "\nUsage: %s <file_list>\n", argv[0]) ;
		fprintf(stderr, "\nOutput: <file_list>-merge.raw\n") ;
		return 1 ;
	}
	
	omp_set_num_threads(16) ;
	
	fsize = 3072 ;
	size = 501 ;
	center = size / 2 ;
	mscale = 1000. ;
	det_flag = 0 ;
	phi = 0. ;
	cenx = 0. ;
	ceny = 0. ;
	detd = 0. ;
	pixsize = 0. ;
	
	model3d = calloc(size*size*size, sizeof(float)) ;
	weight = calloc(size*size*size, sizeof(float)) ;
	frame = malloc(fsize*fsize*sizeof(unsigned short)) ;
	mask = malloc(fsize*fsize*sizeof(int8_t)) ;
	
	// Read in mask
	fp = fopen("data/adsc_mask.byt", "rb") ;
	fread(mask, sizeof(int8_t), fsize*fsize, fp) ;
	fclose(fp) ;
	
	fp_list = fopen(argv[1], "r") ;
	
	while (fscanf(fp_list, "%1023s\n", fname) == 1) {
		fp = fopen(fname, "rb") ;
		
		// Parse header
		fread(header, sizeof(char), 512, fp) ;
		
		char *token = strtok(header, "{};=\n") ;
		x = atoi(strtok(NULL, "{};=\n")) ;
		
		while (token != NULL) {
			if (strcmp(token, "BEAM_CENTER_X") == 0) {
				token = strtok(NULL, "();=\n") ;
				cenx = atof(token) ;
				if (num_files == 0)
					fprintf(stderr, "cenx = %f\n", cenx) ;
			}
			else if (strcmp(token, "BEAM_CENTER_Y") == 0) {
				token = strtok(NULL, "();=\n") ;
				ceny = atof(token) ;
				if (num_files == 0)
					fprintf(stderr, "ceny = %f\n", ceny) ;
			}
			else if (strcmp(token, "PIXEL_SIZE") == 0) {
				token = strtok(NULL, "();=\n") ;
				pixsize = atof(token) ;
				if (num_files == 0)
					fprintf(stderr, "pixsize = %f\n", pixsize) ;
			}
			else if (strcmp(token, "DISTANCE") == 0) {
				token = strtok(NULL, "();=\n") ;
				detd = atof(token) ;
				if (num_files == 0)
					fprintf(stderr, "detd = %f\n", detd) ;
			}
			else if (strcmp(token, "WAVELENGTH") == 0) {
				token = strtok(NULL, "();=\n") ;
				lambda = atof(token) ;
				if (num_files == 0)
					fprintf(stderr, "lambda = %f\n", lambda) ;
			}
			else if (strcmp(token, "PHI") == 0) {
				token = strtok(NULL, "();=\n") ;
				phi = atof(token) ;
				if (num_files == 0)
					fprintf(stderr, "phi = %f\n", phi) ;
				phi *= M_PI / 180. ;
			}
			
			token = strtok(NULL, "{};=\n") ;
		}
		
		if (x > 512)
			fread(header, sizeof(char), x-512, fp) ;
		
		// Parse data
		fread(frame, sizeof(unsigned short), fsize*fsize, fp) ;
		fclose(fp) ;
		
		// For first frame, generate detector
		if (det_flag == 0) {
			cenx /= pixsize ;
			ceny /= pixsize ;
			detd /= pixsize ;
			
			det = malloc(fsize*fsize*3*sizeof(double)) ;
			
			t = 0 ;
			for (x = 0 ; x < fsize ; ++x)
			for (y = 0 ; y < fsize ; ++y) {
				factor = sqrt((x-cenx)*(x-cenx) + (y-ceny)*(y-ceny) + detd*detd) ;
				
				det[t*3 + 0] = mscale * (x-cenx) / factor ;
				det[t*3 + 1] = mscale * (y-ceny) / factor ;
				det[t*3 + 2] = mscale * (detd / factor - 1.) ;
				
				t++ ;
			}
			
			det_flag = 1 ;
			
			fprintf(stderr, "Generated detector\n") ;
		}
		
		// Calculate scale factor
		ave = 0. ;
		for (x = 0 ; x < fsize/3 ; ++x)
		for (y = 0 ; y < fsize/3 ; ++y)
			ave += frame[x*fsize + y] - 40 ;
		for (x = 2*fsize/3 ; x < fsize ; ++x)
		for (y = 0 ; y < fsize/3 ; ++y)
			ave += frame[x*fsize + y] - 40 ;
		for (x = 2*fsize/3 ; x < fsize ; ++x)
		for (y = 2*fsize/3 ; y < fsize ; ++y)
			ave += frame[x*fsize + y] - 40 ;
		for (x = 0 ; x < fsize/3 ; ++x)
		for (y = 2*fsize/3 ; y < fsize ; ++y)
			ave += frame[x*fsize + y] - 40 ;
		
		ave /= 1024*1024*4*20 ;
		
		// Merge in 3D
		make_rot(phi, rot) ;
		
		for (t = 0 ; t < fsize*fsize ; ++t) {
			if (!mask[t])
				continue ;
			
			for (i = 0 ; i < 3 ; ++i) {
				rot_pix[i] = 0. ;
				for (j = 0 ; j < 3 ; ++j)
					rot_pix[i] += rot[i][j] * det[t*3 + j] ;
				rot_pix[i] += center ;
			}
			
			tx = rot_pix[0] ;
			ty = rot_pix[1] ;
			tz = rot_pix[2] ;
			x = tx ;
			y = ty ;
			z = tz ;
			fx = tx - x ;
			fy = ty - y ;
			fz = tz - z ;
			cx = 1. - fx ;
			cy = 1. - fy ;
			cz = 1. - fz ;
			
			w = (frame[t] - 40) / ave ;
			
			f = cx*cy*cz ;
			weight[x*size*size + y*size + z] += f ;
			model3d[x*size*size + y*size + z] += f * w ;
			
			f = cx*cy*fz ;
			weight[x*size*size + y*size + ((z+1)%size)] += f ;
			model3d[x*size*size + y*size + ((z+1)%size)] += f * w ;
			
			f = cx*fy*cz ;
			weight[x*size*size + ((y+1)%size)*size + z] += f ;
			model3d[x*size*size + ((y+1)%size)*size + z] += f * w ;
			
			f = cx*fy*fz ;
			weight[x*size*size + ((y+1)%size)*size + ((z+1)%size)] += f ;
			model3d[x*size*size + ((y+1)%size)*size + ((z+1)%size)] += f * w ;
			
			f = fx*cy*cz ;
			weight[((x+1)%size)*size*size + y*size + z] += f ;
			model3d[((x+1)%size)*size*size + y*size + z] += f * w ;
			
			f = fx*cy*fz ;
			weight[((x+1)%size)*size*size + y*size + ((z+1)%size)] += f ;
			model3d[((x+1)%size)*size*size + y*size + ((z+1)%size)] += f * w ;
			
			f = fx*fy*cz ;
			weight[((x+1)%size)*size*size + ((y+1)%size)*size + z] += f ;
			model3d[((x+1)%size)*size*size + ((y+1)%size)*size + z] += f * w ;
			
			f = fx*fy*fz ;
			weight[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += f ;
			model3d[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += f * w ;
		}
		
		num_files++ ;
		fprintf(stderr, "\rMerged %.4d files", num_files) ;
	}
	fprintf(stderr, "\n") ;
	
	fclose(fp_list) ;
	
	// Normalize by interpolation weights
	for (t = 0 ; t < size*size*size ; ++t)
	if (weight[t] != 0.)
		model3d[t] /= weight[t] ;
	
	// Symmetrize merge
	for (x = 0 ; x <= center ; ++x)
	for (y = 0 ; y <= center ; ++y)
	for (z = 0 ; z <= center ; ++z) {
		ave = 0. ;
		avecount = 0 ;
		
		vox = x*size*size + y*size + z ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = x*size*size + y*size + (2*center-z) ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = x*size*size + (2*center-y)*size + z ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = x*size*size + (2*center-y)*size + (2*center-z) ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = (2*center-x)*size*size + y*size + z ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = (2*center-x)*size*size + y*size + (2*center-z) ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = (2*center-x)*size*size + (2*center-y)*size + z ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		vox = (2*center-x)*size*size + (2*center-y)*size + (2*center-z) ;
		if (model3d[vox] > 0.) {
			ave += model3d[vox] ;
			avecount++ ;
		}
		
		if (avecount > 0)
			ave /= avecount ;
		
		model3d[x*size*size + y*size + z] = ave ;
		model3d[x*size*size + y*size + (2*center-z)] = ave ;
		model3d[x*size*size + (2*center-y)*size + z] = ave ;
		model3d[x*size*size + (2*center-y)*size + (2*center-z)] = ave ;
		model3d[(2*center-x)*size*size + y*size + z] = ave ;
		model3d[(2*center-x)*size*size + y*size + (2*center-z)] = ave ;
		model3d[(2*center-x)*size*size + (2*center-y)*size + z] = ave ;
		model3d[(2*center-x)*size*size + (2*center-y)*size + (2*center-z)] = ave ;
	}
	
	// Write to file
	sprintf(fname, "%s_merge.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "Writing merge to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model3d, sizeof(float), size*size*size, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(mask) ;
	free(det) ;
	free(model3d) ;
	free(weight) ;
	free(frame) ;
	
	return 0 ;
}
