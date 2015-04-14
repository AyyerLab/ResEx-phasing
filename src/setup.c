#include "brcont.h"

int allocate_memory(int) ;
int parse_intens(char*) ;
int parse_hkl(char*, double) ;
int parse_support(char*) ;
void create_plans(char*) ;
int gen_input(char*, int) ;
int parse_quat(char*) ;

int setup() {
	char var_name[500], input_fname[500], quat_fname[500] ;
	char intens_fname[500], hkl_fname[500], support_fname[500], wisdom_fname[500] ;
	double hklqmax ;
	
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
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", input_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", quat_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", wisdom_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%lf\n", &hklqmax) ;
	fclose(fp) ;
	
	vol = size*size*size ;
	hklvol = hsize*ksize*lsize ;
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(4) ;
	
	if (allocate_memory(1))
		return 1 ;
	if (parse_intens(intens_fname))
		return 1 ;
	if (parse_hkl(hkl_fname, hklqmax))
		return 1 ;
	if (parse_support(support_fname))
		return 1 ;
	if (parse_quat(quat_fname))
		return 1 ;
	gen_input(input_fname, 0) ;
	create_plans(wisdom_fname) ;
	
	return 0 ;
}	

int setup_gen() {
	char var_name[500], input_fname[500], quat_fname[500] ;
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
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", input_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", quat_fname) ;
	fgets(var_name, 500, fp) ; fscanf(fp, "%s\n", wisdom_fname) ;
	fclose(fp) ;
	
	vol = size*size*size ;
	hklvol = hsize*ksize*lsize ;
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(4) ;
	
	if (allocate_memory(0))
		return 1 ;
	if (gen_input(input_fname, 1))
		return 1 ;
	if (parse_quat(quat_fname))
		return 1 ;
	create_plans(wisdom_fname) ;
	
	return 0 ;
}

int allocate_memory(int flag) {
	iterate = malloc(vol * sizeof(float)) ;
	exp_mag = malloc(vol * sizeof(float)) ;
	
	if (flag == 1) {
		obs_mag = malloc(vol * sizeof(float)) ;
		hkl_calc = fftwf_malloc(hklvol * sizeof(fftwf_complex)) ;
		p1 = malloc(vol * sizeof(float)) ;
		p2 = malloc(vol * sizeof(float)) ;
		r1 = malloc(vol * sizeof(float)) ;
	}
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	support = malloc(vol * sizeof(long)) ;
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
				= sqrt(intens[i*s*s + j*s + k]) * 370.67506 ;
		else
			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= intens[i*s*s + j*s + k] ;
	
	free(intens) ;
	
	return 0 ;
}

int parse_hkl(char *fname, double hklqmax) {
	long h, k, l ;
	long hc = hsize/2, kc = ksize/2, lc = lsize/2 ;
	double hc2 = hc*hc, kc2 = kc*kc, lc2 = lc*lc ;
	double dist ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found. Exiting.\n", fname) ;
		return 1 ;
	}
	
	fftwf_complex *hkl_temp = malloc(hklvol * sizeof(fftw_complex)) ;
	
	fread(hkl_temp, sizeof(fftw_complex), hklvol, fp) ;
	fclose(fp) ;
	
	for (h = 0 ; h < hsize ; ++h)
	for (k = 0 ; k < ksize ; ++k)
	for (l = 0 ; l < lsize ; ++l) {
		dist = (h-hc)*(h-hc)/hc2 + (k-kc)*(k-kc)/kc2 + (l-lc)*(l-lc)/lc2 ;
		
		if (dist < hklqmax*hklqmax)
			hkl_calc[((h+hc+1)%hsize)*ksize*lsize + ((k+kc+1)%ksize)*lsize + ((l+lc+1)%lsize)] 
			 = hkl_temp[h*ksize*lsize + k*lsize + l] 
			   * pow(-1., h-hc+k-kc+l-lc)
			   * cexp(-4*I*M_PI*(l-lc)/lsize) ;
		else
			hkl_calc[((h+hc+1)%hsize)*ksize*lsize + ((k+kc+1)%ksize)*lsize + ((l+lc+1)%lsize)] = FLT_MAX ;
	}
	
	fftwf_free(hkl_temp) ;
	
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
	
	forward_hkl = fftwf_plan_dft_3d(hsize, ksize, lsize, rhkl, fhkl, FFTW_FORWARD, FFTW_MEASURE) ;
	inverse_hkl = fftwf_plan_dft_3d(hsize, ksize, lsize, fhkl, rhkl, FFTW_BACKWARD, FFTW_MEASURE) ;
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

int parse_quat(char *fname) {
	long r, t ;
	double tot_weight = 0. ;
	
	
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found.\n", fname) ;
		return 1 ;
	}
	fscanf(fp, "%d", &num_rot) ;
	quat = malloc(num_rot * 5 * sizeof(double)) ;
	for (r = 0 ; r < num_rot ; ++r) {
		for (t = 0 ; t < 5 ; ++t)
			fscanf(fp, "%lf ", &quat[r*5 + t]) ;
		tot_weight += quat[r*5 + 4] ;
	}
	
	for (r = 0 ; r < num_rot ; ++r)
		quat[r*5 + 4] /= tot_weight ;
	
	fclose(fp) ;
	
	return 0 ;
}
