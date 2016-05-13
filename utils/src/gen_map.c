#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	long s, vol, mvol, x, y, z, vox, num_supp = 0 ;
	long mmin_x, mmax_x, ms_x ;
	long mmin_y, mmax_y, ms_y ;
	long mmin_z, mmax_z, ms_z ;
	int ival ;
	float fval, vox_x, vox_y, vox_z, *model, *mmodel ;
	float mean = 0.f, rms = 0.f, max = -1.e20, min = 1.e20 ;
	char *buffer, fname[999] ;
	uint8_t *support ;
	FILE *fp ;
	
	if (argc < 7) {
		fprintf(stderr, "Format: %s <recon_fname> <size> <vox_x> <vox_y> <vox_z> <supp_fname>\n", argv[0]) ;
		return 1 ;
	}
	s = atoi(argv[2]) ;
	vox_x = atof(argv[3]) ;
	vox_y = atof(argv[4]) ;
	vox_z = atof(argv[5]) ;
	vol = s*s*s ;
	
	// Read model and support
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	support = malloc(vol * sizeof(uint8_t)) ;
	fp = fopen(argv[6], "rb") ;
	fread(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;

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
		mmodel[((z+ms_z/2)%ms_z)*ms_y*ms_x + ((y+ms_y/2)%ms_y)*ms_x + ((x+ms_x/2)%ms_x)] // ms_x*ms_y*ms_z, translated
//		mmodel[z*ms_y*ms_x + y*ms_x + x] // ms_x*ms_y*ms_z
			= model[(x+mmin_x)*s*s + (y+mmin_y)*s + (z+mmin_z)] ; // s*s*s
	free(model) ;

	// Start writing map
	// --------------------------------------------------------------------------------
	sprintf(fname, "data/maps/%s.map.ccp4", remove_ext(extract_fname(argv[1]))) ;
	fp = fopen(fname, "wb") ;
	// NC, NR, NS
	ival = ms_x ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_y ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_z ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// MODE
	ival = 2 ;
	fwrite(&ival, sizeof(int), 1, fp) ;

	// NCSTART, NRSTART, NSSTART
	ival = 0 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;

	// NX, NY, NZ
	ival = ms_z ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_y ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_x ;
	fwrite(&ival, sizeof(int), 1, fp) ;

	// X_LENGTH, Y_LENGTH, Z_LENGTH
	fval = vox_z * ms_z ;
	fprintf(stderr, "box size = %f", fval) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fval = vox_y * ms_y ;
	fprintf(stderr, " x %f", fval) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fval = vox_x * ms_x ;
	fprintf(stderr, " x %f A\n", fval) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	// ALPHA, BETA, GAMMA
	fval = 90. ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fwrite(&fval, sizeof(float), 1, fp) ;

	// MAPC, MAPR, MAPS
	ival = 3 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = 2 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = 1 ;
	fwrite(&ival, sizeof(int), 1, fp) ;

	// AMIN, AMAX, AMEAN
	fwrite(&max, sizeof(float), 1, fp) ;
	fwrite(&min, sizeof(float), 1, fp) ;
	fwrite(&mean, sizeof(float), 1, fp) ;

	// ISPG, NSYMBT
	ival = 1 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = 0 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// LSKFlG, SKWMAT, SKWTRN
	ival = 0 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fval = 0. ;
	for (x = 0 ; x < 12 ; ++x)
		fwrite(&fval, sizeof(float), 1, fp) ;

	// EXTRA
	buffer = calloc(15*4, sizeof(char)) ;
	fwrite(buffer, sizeof(char), 15*4, fp) ;
	free(buffer) ;
	// MAP
	buffer = malloc(4 * sizeof(char)) ;
	buffer[0] = 'M' ;
	buffer[1] = 'A' ;
	buffer[2] = 'P' ;
	buffer[3] = ' ' ;
	fwrite(buffer, sizeof(char), 4, fp) ;
	free(buffer) ;

	// MACHST
	buffer = malloc(4 * sizeof(char)) ;
	buffer[0] = 0x44 ;
	buffer[1] = 0x41 ;
	buffer[2] = 0x00 ;
	buffer[3] = 0x00 ;
	fwrite(buffer, sizeof(char), 4, fp) ;
	free(buffer) ;

	// RMS
	fwrite(&rms, sizeof(float), 1, fp) ;
	// NLABL
	ival = 1 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// LABEL_N
	char label[800] ;
	sprintf(label, "%s %ld %s", argv[1], s, argv[6]) ;
	fprintf(stderr, "label = %s\n", label) ;
	fwrite(label, sizeof(char), 800, fp) ;

	// VOXELS
	fwrite(mmodel, sizeof(float), mvol, fp) ;
	
	fclose(fp) ;
	
	free(mmodel) ;
	
	return 0 ;
}
