#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
	long s, vol, ms, mvol, x, y, z ;
	int ival ;
	float fval, *model, *mmodel ;
	float mean = 0.f, rms = 0.f, max = -1.e20, min = 1.e20 ;
	char *buffer, fname[999] ;
	FILE *fp ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <recon_fname>\n", argv[0]) ;
		return 1 ;
	}
	
	s = 501 ;
	vol = s*s*s ;
	ms = 171 ;
	mvol = ms*ms*ms ;
	
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	mmodel = malloc(mvol * sizeof(float)) ;
	for (x = 163 ; x < 334 ; ++x)
	for (y = 163 ; y < 334 ; ++y)
	for (z = 163 ; z < 334 ; ++z) {
		fval = model[x*s*s + y*s + z] ;
//		if (fval < 0.5)
//			fval = 0. ;
//		if (fval > 10.)
//			fval = 10. ;
		mmodel[(z-163)*ms*ms + (y-163)*ms + (x-163)] = fval ;
	}
	free(model) ;
	
	for (x = 0 ; x < mvol ; ++x) {
		mean += mmodel[x] ;
		rms += mmodel[x]*mmodel[x] ;
		if (mmodel[x] > max)
			max = mmodel[x] ;
		if (mmodel[x] < min)
			min = mmodel[x] ;
	}
	mean /= (float) mvol ;
	rms /= (float) mvol ;
	rms = sqrtf(rms - mean*mean) ;
	
	fp = fopen("data/mmodel.map", "wb") ;
	fwrite(mmodel, sizeof(float), mvol, fp) ;
	fclose(fp) ;
	
	fprintf(stderr, "max = %.6e, min = %.6e\n", max, min) ;
	fprintf(stderr, "mean = %.6e, rms = %.6e\n", mean, rms) ;
	
	sprintf(fname, "data/maps/%s.map.ccp4", remove_ext(extract_fname(argv[1]))) ;
	fp = fopen(fname, "wb") ;
	// NC, NR, NS
	ival = ms ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;
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
	ival = ms ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	fwrite(&ival, sizeof(int), 1, fp) ;
	// X_LENGTH, Y_LENGTH, Z_LENGTH
//	fval = 542.07 ;
//	fval = 271.035 ;
//	fval = 272.883 ;
//	fval = 271.83 ;
	fval = 315.07 ;
	fwrite(&fval, sizeof(float), 1, fp) ;
//	fval = 624.5 ;
//	fval = 312.25 ;
//	fval = 308.906 ;
//	fval = 307.71 ;
	fval = 309.52 ;
	fwrite(&fval, sizeof(float), 1, fp) ;
//	fval = 632.7 ;
//	fval = 316.35 ;
//	fval = 314.445 ;
//	fval = 313.23 ;
	fval = 273.49 ;
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
