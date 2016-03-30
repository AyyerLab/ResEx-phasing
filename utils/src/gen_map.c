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
	long s, vol, mvol, x, y, z, num_supp ;
	long mmin_x, mmax_x, ms_x ;
	long mmin_y, mmax_y, ms_y ;
	long mmin_z, mmax_z, ms_z ;
	int ival ;
	float fval, voxres, *model, *mmodel ;
	float mean = 0.f, rms = 0.f, max = -1.e20, min = 1.e20 ;
	char *buffer, fname[999] ;
	uint8_t *support, *msupp ;
	FILE *fp ;
	
	if (argc < 5) {
		fprintf(stderr, "Format: %s <recon_fname> <size> <voxres> <supp_fname>\n", argv[0]) ;
		return 1 ;
	}
	s = atoi(argv[2]) ;
	voxres = atof(argv[3]) ;
	vol = s*s*s ;
	
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	support = malloc(vol * sizeof(uint8_t)) ;
	fp = fopen(argv[4], "rb") ;
	fread(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	mmin_x = s ; mmin_y = s ; mmin_z = s ;
	mmax_x = 0 ; mmax_y = 0 ; mmax_z = 0 ;
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z) {
		if (support[x*s*s + y*s + z] == 0)
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
	}
	mmin_x -= 2 ; mmin_y -= 2 ; mmin_z -= 2 ;
	mmax_x += 2 ; mmax_y += 2 ; mmax_z += 2 ;
//	mmin_x = 127, mmin_y = 127, mmin_z = 127 ;
//	mmax_x = 175, mmax_y = 175, mmax_z = 175 ;
	ms_x = mmax_x - mmin_x + 1 ;
	ms_y = mmax_y - mmin_y + 1 ;
	ms_z = mmax_z - mmin_z + 1 ;
	mvol = ms_x * ms_y * ms_z ;
	fprintf(stderr, "min = (%ld, %ld, %ld)\n", mmin_x, mmin_y, mmin_z) ;
	fprintf(stderr, "max = (%ld, %ld, %ld)\n", mmax_x, mmax_y, mmax_z) ;
	fprintf(stderr, "Model volume = %ld x %ld x %ld = %ld\n", ms_x, ms_y, ms_z, mvol) ;
		
	mmodel = malloc(mvol * sizeof(float)) ;
	msupp = malloc(mvol * sizeof(uint8_t)) ;
	for (x = mmin_x ; x <= mmax_x ; ++x)
	for (y = mmin_y ; y <= mmax_y ; ++y)
	for (z = mmin_z ; z <= mmax_z ; ++z) {
		mmodel[(x-mmin_x)*ms_y*ms_z + (y-mmin_y)*ms_z + (z-mmin_z)] = model[x*s*s + y*s + z] ;
		msupp[(x-mmin_x)*ms_y*ms_z + (y-mmin_y)*ms_z + (z-mmin_z)] = support[x*s*s + y*s + z] ;
	}
	free(model) ;
	free(support) ;
	
	num_supp = 0 ;
	for (x = 0 ; x < mvol ; ++x) {
		if (msupp[x] == 0)
			continue ;
		
		mean += mmodel[x] ;
		rms += mmodel[x]*mmodel[x] ;
		if (mmodel[x] > max)
			max = mmodel[x] ;
		if (mmodel[x] < min)
			min = mmodel[x] ;
		
		num_supp++ ;
	}
	mean /= (float) num_supp ;
	rms /= (float) num_supp ;
	rms = sqrtf(rms - mean*mean) ;
	
	fprintf(stderr, "Effective solvent fraction = %f\n", (double) num_supp / mvol) ;
	fprintf(stderr, "max = %.6e, min = %.6e\n", max, min) ;
	fprintf(stderr, "mean = %.6e, rms = %.6e\n", mean, rms) ;
	
	sprintf(fname, "data/maps/%s.map.ccp4", remove_ext(extract_fname(argv[1]))) ;
	fp = fopen(fname, "wb") ;
	// NC, NR, NS
	ival = ms_z ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_y ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_x ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// MODE
	ival = 2 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// NCSTART, NRSTART, NSSTART
//	ival = -85 ;
	ival = 0 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// NX, NY, NZ
	ival = ms_x ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_y ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = ms_z ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// X_LENGTH, Y_LENGTH, Z_LENGTH
	fval = voxres/(s-1) * ms_x ;
	fprintf(stderr, "box size = %f", fval) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fval = voxres/(s-1) * ms_y ;
	fprintf(stderr, " x %f", fval) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fval = voxres/(s-1) * ms_z ;
	fprintf(stderr, " x %f A\n", fval) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	// ALPHA, BETA, GAMMA
	fval = 90. ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	fwrite(&fval, sizeof(float), 1, fp) ;
	// MAPC, MAPR, MAPS
	ival = 1 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = 2 ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	ival = 3 ;
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
	char label[800] = "::CFEL DataBank::0001::" ;
	fwrite(label, sizeof(char), 800, fp) ;
	
	// VOXELS
	fwrite(mmodel, sizeof(float), mvol, fp) ;
	
	fclose(fp) ;
	
	free(mmodel) ;
	
	return 0 ;
}
