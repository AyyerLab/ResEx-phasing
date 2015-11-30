#include "brcont.h"

int allocate_memory(int) ;
int parse_intens(char*) ;
int parse_bragg(char*, double) ;
int parse_support(char*) ;
void create_plans(char*) ;
int gen_input(char*) ;
int parse_quat(char*) ;
void gen_mask(float) ;

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
	fftwf_plan_with_nthreads(32) ;
	
	if (allocate_memory(1))
		return 1 ;
	if (parse_intens(intens_fname))
		return 1 ;
	if (parse_support(support_fname))
		return 1 ;
	gen_mask(300.f) ;
	gen_input(input_fname) ;
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
	fftwf_plan_with_nthreads(32) ;
	
	if (allocate_memory(0))
		return 1 ;
	if (parse_bragg(bragg_fname, 1.))
		return 1 ;
	
	return 0 ;
}

int allocate_memory(int flag) {
	incoh_mag = malloc(vol * sizeof(float)) ;
	coh_mag = malloc(vol * sizeof(float)) ;
	
	if (flag == 0) {
		intens = calloc(vol, sizeof(float)) ;
		bragg_calc = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	}
	if (flag == 1) {
		p1 = malloc(vol * sizeof(float)) ;
		p2 = malloc(vol * sizeof(float)) ;
		r1 = malloc(vol * sizeof(float)) ;
		iterate = malloc(vol * sizeof(float)) ;
		support = malloc(vol * sizeof(long)) ;
		obs_mag = malloc(vol * sizeof(float)) ;
		mask = calloc(vol, sizeof(uint8_t)) ;
		rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
		fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	}
	
	return 0 ;
}

int parse_intens(char *fname) {
	long i, j, k, s = size ;//, c = size/2 ;
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
//			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
			obs_mag[i*s*s + j*s + k]
				= sqrt(intens[i*s*s + j*s + k]) ;
		else
//			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
			obs_mag[i*s*s + j*s + k]
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
//			 = bragg_temp[x*size*size + y*size + z] 
//			   * powf(-1.f, x-c+y-c+z-c) ;
			 = bragg_temp[x*size*size + y*size + z] ;
		else
//			bragg_calc[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = FLT_MAX ;
			bragg_calc[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = 0. ;
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

int gen_input(char *fname) {
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Random start\n") ;
		init_model(iterate) ;
	}
	else {
		fprintf(stderr, "Starting from %s\n", fname) ;
		fread(iterate, sizeof(float), vol, fp) ;
		fclose(fp) ;
	}
	
	return 0 ;
}

void gen_mask(float rmax) {
	long x, y, z ;
	long s = size, c = s/2 ;
	double dist ;
	
	// Default value for mask is NON_BRAGG
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z) {
		dist = sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		
		if (dist > c) // Ignore cube corners
			mask[x*s*s + y*s + z] = IGNORE ;
		else if (dist > rmax) // High resolution is all non-Bragg
			mask[x*s*s + y*s + z] = NON_BRAGG ;
//		else if ((x-c)%3 == 0 && (y-c)%4 == 0 && (z-c)%7 == 0) {
//			mask[x*s*s + y*s + z] = BRAGG ;
/*			mask[(x+1)*s*s + y*s + z] = IGNORE ;
			mask[x*s*s + (y+1)*s + z] = IGNORE ;
			mask[x*s*s + y*s + (z+1)] = IGNORE ;
			mask[(x-1)*s*s + y*s + z] = IGNORE ;
			mask[x*s*s + (y-1)*s + z] = IGNORE ;
			mask[x*s*s + y*s + (z-1)] = IGNORE ;
*///		}
	}
	
	FILE *fp = fopen("data/mask.supp", "wb") ;
	fwrite(mask, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
}
