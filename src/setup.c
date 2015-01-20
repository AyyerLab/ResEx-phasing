#include "brcont.h"

int allocate_memory() ;
int parse_intens(char*) ;
int parse_hkl(char*) ;
int parse_support(char*) ;
void create_plans(char*) ;

int setup() {
	char var_name[500] ;
	char intens_fname[500], hkl_fname[500], support_fname[500], wisdom_fname[500] ;
	
	FILE *fp = fopen("src/config.conf", "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Config file not found\n") ;
		return 1 ;
	}
	fgets(var_name, 500, fp) ; fscanf(fp, "%ld\n", &size) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%ld %ld %ld\n", &hsize, &ksize, &lsize) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%ld %ld %ld\n", &hoffset, &koffset, &loffset) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", intens_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", hkl_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", support_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", wisdom_fname) ;
	fclose(fp) ;
	
	vol = size*size*size ;
	hklvol = hsize*ksize*lsize ;
	
	if (allocate_memory())
		return 1 ;
	if (parse_intens(intens_fname))
		return 1 ;
	if (parse_hkl(hkl_fname))
		return 1 ;
	if (parse_support(support_fname))
		return 1 ;
	create_plans(wisdom_fname) ;
	
	return 0 ;
}	

int allocate_memory() {
	int i ;
	
	obs_mag = malloc(vol * sizeof(float)) ;
	exp_intens = malloc(vol * sizeof(float)) ;
	hkl_mag = malloc(hklvol * sizeof(float)) ;
	exp_hkl = malloc(hklvol * sizeof(float)) ;
	
	for (i = 0 ; i < 3 ; ++i) {
		p1[i] = malloc(vol * sizeof(float)) ;
		p2[i] = malloc(vol * sizeof(float)) ;
		r1[i] = malloc(vol * sizeof(float)) ;
	}
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	rhkl = fftwf_malloc(hklvol * sizeof(fftwf_complex)) ;
	fhkl = fftwf_malloc(hklvol * sizeof(fftwf_complex)) ;
	
	return 0 ;
}

int parse_intens(char *fname) {
	long i, j, k, s = size, c = size/2 ;
	float *intens ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found. Exiting.\n", fname) ;
		return 1 ;
	}
	
	intens = malloc(vol * sizeof(float)) ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	for (i = 0 ; i < s ; ++i)
	for (j = 0 ; j < s ; ++j)
	for (k = 0 ; k < s ; ++k)
		if (intens[i*s*s + j*s + k] > 0.)
			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= sqrt(intens[i*s*s + j*s + k]) ;
		else
			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= intens[i*s*s + j*s + k] ;
	
	free(intens) ;
	
	return 0 ;
}

int parse_hkl(char *fname) {
	long i, j, k ;
	long hs = hsize, hc = hsize/2 ;
	long ks = ksize, kc = ksize/2 ;
	long ls = lsize, lc = lsize/2 ;
	float *intens ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found. Exiting.\n", fname) ;
		return 1 ;
	}
	
	intens = malloc(vol * sizeof(float)) ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	for (i = 0 ; i < hs ; ++i)
	for (j = 0 ; j < ks ; ++j)
	for (k = 0 ; k < ls ; ++k)
		if (intens[i*ks*ls + j*ls + k] > 0.)
			hkl_mag[((i+hc+1)%hs)*ks*ls + ((j+kc+1)%hs)*ls + ((k+lc+1)%hs)]
				= sqrt(intens[i*ks*ls + j*ls + k]) ;
		else
			hkl_mag[((i+hc+1)%hs)*ks*ls + ((j+kc+1)%ks)*ls + ((k+lc+1)%ls)]
				= intens[i*ks*ls + j*ls + k] ;
	
	free(intens) ;
	
	return 0 ;
}

int parse_support(char *fname) {
	long x, t ;
	uint8_t *supvol = malloc(vol * sizeof(uint8_t)) ;
	
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found.\n", fname) ;
		return 1 ;
	}
	fread(supvol, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	num_supp = 0 ;
	for (x = 0 ; x < vol ; ++x)
		num_supp += supvol[x] ;
	
	support = malloc(num_supp * sizeof(long)) ;
	t = 0 ;
	for (x = 0 ; x < vol ; ++x)
	if (supvol[x] == 1)
		support[t++] = x ;
	
	free(supvol) ;
	
	return 0 ;
}

void create_plans(char *fname) {
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Measuring plans\n") ;
		forward_cont = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse_cont = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
		
		fp = fopen(fname, "wb") ;
		fftwf_export_wisdom_to_file(fp) ;
		fclose(fp) ;
	
		fprintf(stderr, "Created plans\n") ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward_cont = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse_cont = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}
	
	fftwf_plan_dft_3d(hsize, ksize, lsize, rhkl, fhkl, FFTW_FORWARD, FFTW_MEASURE) ;
	fftwf_plan_dft_3d(hsize, ksize, lsize, fhkl, rhkl, FFTW_BACKWARD, FFTW_MEASURE) ;
}
