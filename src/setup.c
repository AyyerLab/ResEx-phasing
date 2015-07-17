#include "brcont.h"

int allocate_memory(int) ;
int parse_intens(char*) ;
int parse_bragg(char*, double) ;
int parse_support(char*) ;
void create_plans(char*) ;
int gen_input(char*, int) ;
int parse_quat(char*) ;

int setup() {
	char var_name[500], input_fname[500] ;
	char intens_fname[500], bragg_fname[500], support_fname[500], wisdom_fname[500] ;
	double braggqmax ;
	
	FILE *fp = fopen("src/config.conf", "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Config file not found\n") ;
		return 1 ;
	}
	fgets(var_name, 500, fp) ; fscanf(fp, "%ld\n", &size) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", intens_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", bragg_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", support_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", input_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", wisdom_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%lf\n", &braggqmax) ;
	fclose(fp) ;
	
	vol = size*size*size ;
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	if (allocate_memory(1))
		return 1 ;
	if (parse_intens(intens_fname))
		return 1 ;
	if (parse_bragg(bragg_fname, braggqmax))
		return 1 ;
	if (parse_support(support_fname))
		return 1 ;
	gen_input(input_fname, 0) ;
	create_plans(wisdom_fname) ;
	
	return 0 ;
}	

int setup_gen() {
	char var_name[500], input_fname[500] ;
	char intens_fname[500], bragg_fname[500], support_fname[500], wisdom_fname[500] ;
	
	FILE *fp = fopen("src/config.conf", "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Config file not found\n") ;
		return 1 ;
	}
	fgets(var_name, 500, fp) ; fscanf(fp, "%ld\n", &size) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", intens_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", bragg_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", support_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", input_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", wisdom_fname) ;
	fclose(fp) ;
	
	vol = size*size*size ;
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(16) ;
	
	if (allocate_memory(0))
		return 1 ;
	if (gen_input(input_fname, 1))
		return 1 ;
	create_plans(wisdom_fname) ;
	
	return 0 ;
}

int allocate_memory(int flag) {
	iterate = malloc(vol * sizeof(float)) ;
	exp_mag = malloc(vol * sizeof(float)) ;
	
	if (flag == 1) {
		obs_mag = malloc(vol * sizeof(float)) ;
		bragg_calc = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
		p1 = malloc(vol * sizeof(float)) ;
		p2 = malloc(vol * sizeof(float)) ;
		r1 = malloc(vol * sizeof(float)) ;
	}
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	support = malloc(vol * sizeof(long)) ;
	
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
//				= sqrt(intens[i*s*s + j*s + k]) * 370.67506 ; // anton_50pr_subt
//				= sqrt(intens[i*s*s + j*s + k]) * 157.19 ; // lorenzo_2848_hard
//				= sqrt(intens[i*s*s + j*s + k]) * 421.17 ; // lorenzo_hard_iso_sym
//				= sqrt(intens[i*s*s + j*s + k]) * 343.06 ; // lorenzo_hard_iso_sym w/ dominik-3
//				= sqrt(intens[i*s*s + j*s + k]) * 341.33 ; // lorenzo_hard_iso_sym w/ dominik-3-test
				= sqrt(intens[i*s*s + j*s + k]) * 336.16 ; // lorenzo_hard_iso_sym w/ dominik-4
		else
			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= intens[i*s*s + j*s + k] ;
	
	free(intens) ;
	
	return 0 ;
}

int parse_bragg(char *fname, double braggqmax) {
	long x, y, z, c = size/2 ;
	double dist, c2 = c*c ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found. Exiting.\n", fname) ;
		return 1 ;
	}
	
	fftwf_complex *bragg_temp = malloc(vol * sizeof(fftw_complex)) ;
	
	fread(bragg_temp, sizeof(fftw_complex), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dist = ((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) / c2 ;
		
		if (dist < braggqmax*braggqmax) 
			bragg_calc[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] 
//			 = bragg_temp[z*size*size + y*size + x] // Reversing axes
			 = bragg_temp[x*size*size + y*size + z] 
			   * powf(-1.f, x-c+y-c+z-c) ;
		else
			bragg_calc[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = FLT_MAX ;
	}
	
	fftwf_free(bragg_temp) ;
	
	return 0 ;
}

int parse_support(char *fname) {
	long x ;
	uint8_t *supvol = malloc(vol * sizeof(uint8_t)) ;
	
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found.\n", fname) ;
		return 1 ;
	}
	fread(supvol, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < vol ; ++x)
	if (supvol[x] == 1)
		support[num_supp++] = x ;
	
	free(supvol) ;
	
	return 0 ;
}

void create_plans(char *fname) {
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Measuring plans\n") ;
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
		
		fp = fopen(fname, "wb") ;
		fftwf_export_wisdom_to_file(fp) ;
		fclose(fp) ;
	
		fprintf(stderr, "Created plans\n") ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}
}

int gen_input(char *fname, int flag) {
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		if (flag == 0) {
			fprintf(stderr, "Random start\n") ;
			init_model(iterate) ;
		}
		else {
			fprintf(stderr, "Cannot find input %s\n", fname) ;
			return 1 ;
		}
	}
	else {
		fprintf(stderr, "Starting from %s\n", fname) ;
		fread(iterate, sizeof(float), vol, fp) ;
		fclose(fp) ;
	}
	
	return 0 ;
}
