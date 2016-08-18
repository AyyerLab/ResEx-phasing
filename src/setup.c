#include "brcont.h"

int allocate_memory(int) ;
int parse_intens(char*, float) ;
int parse_bragg(char*, double) ;
int parse_support(char*) ;
void create_plans(char*) ;
int gen_input(char*, int) ;
int parse_quat(char*) ;
int read_histogram(char*, long) ;

int setup(char *config_fname) {
	char line[999], *token ;
	char input_fname[999], hist_fname[999] ;
	char intens_fname[999], bragg_fname[999] ;
	char support_fname[999], wisdom_fname[999] ;
	double bragg_qmax = 0. ;
	float scale_factor = 0. ;
	int num_threads = -1 ;
	
	size = 0 ;
	output_prefix[0] = '\0' ;
	strcpy(point_group, "222") ;
	do_histogram = 0 ;
	do_positivity = 0 ;
	strcpy(algorithm_name, "DM") ;
	strcpy(avg_algorithm_name, "Avg") ;
	
	FILE *fp = fopen(config_fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Config file, %s, not found.\n", config_fname) ;
		return 1 ;
	}
	while (fgets(line, 999, fp) != NULL) {
		token = strtok(line, " =") ;
		if (token[0] == '#' || token[0] == '\n' || token[0] == '[')
			continue ;
		
		// Parameters
		if (strcmp(token, "size") == 0)
			size = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "bragg_qmax") == 0)
			bragg_qmax = atof(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "scale_factor") == 0)
			scale_factor = atof(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "num_threads") == 0)
			num_threads = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "point_group") == 0)
			strcpy(point_group, strtok(NULL, " =\n")) ;
		// Files
		else if (strcmp(token, "intens_fname") == 0)
			strcpy(intens_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "bragg_fname") == 0)
			strcpy(bragg_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "input_fname") == 0)
			strcpy(input_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "support_fname") == 0)
			strcpy(support_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "output_prefix") == 0)
			strcpy(output_prefix, strtok(NULL, " =\n")) ;
		// Algorithm
		else if (strcmp(token, "algorithm") == 0)
			strcpy(algorithm_name, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "avg_algorithm") == 0)
			strcpy(avg_algorithm_name, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "beta") == 0)
			algorithm_beta = atof(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "histogram") == 0)
			do_histogram = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "positivity") == 0)
			do_positivity = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "hist_fname") == 0)
			strcpy(hist_fname, strtok(NULL, " =\n")) ;
	}
	
	if (size == 0) {
		fprintf(stderr, "Need nonzero size in config file\n") ;
		return 1 ;
	}
	else if (size%2 == 0)
		fprintf(stderr, "size = %ld is even. An odd number is preferred.\n", size) ;
	
	if (strcmp(point_group, "222")*strcmp(point_group, "4")*strcmp(point_group, "1") != 0) {
		fprintf(stderr, "Point group needs to be either 1, 222 or 4 for now\n") ;
		return 1 ;
	}
	
	if (output_prefix[0] == '\0') {
		fprintf(stderr, "Using default output prefix data/output\n") ;
		strcpy(output_prefix, "data/output") ;
	}
	
	if (num_threads == -1)
		num_threads = omp_get_max_threads() ;
	sprintf(wisdom_fname, "data/wisdom_%ld_%d", size, num_threads) ;
	
	fprintf(stderr, "Scale factor = %f\n", scale_factor) ;
	fprintf(stderr, "Symmetrizing with point group: %s\n", point_group) ;
	
	vol = size*size*size ;
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(num_threads) ;
	
	if (allocate_memory(1))
		return 1 ;
	if (parse_intens(intens_fname, scale_factor))
		return 1 ;
	if (parse_bragg(bragg_fname, bragg_qmax))
		return 1 ;
	if (parse_support(support_fname))
		return 1 ;
	if (do_histogram) {
		if (read_histogram(hist_fname, num_supp))
			return 1 ;
	}
	gen_input(input_fname, 0) ;
	create_plans(wisdom_fname) ;
	
	sprintf(line, "%s-log.dat", output_prefix) ;
	fp = fopen(line, "w") ;
	fprintf(fp, "Resolution extension iterative phasing\n") ;
	fprintf(fp, "Data: %s %s\n", bragg_fname, intens_fname) ;
	fprintf(fp, "Support: %s (%ld)\n", support_fname, num_supp) ;
	fprintf(fp, "Algorithm: %s with beta = %.2f\n", algorithm_name, algorithm_beta) ;
	fprintf(fp, "Averaging algorithm: %s\n", avg_algorithm_name) ;
	if (do_histogram)
		fprintf(fp, "Applying histogram constraint: %s\n", hist_fname) ;
	else
		fprintf(fp, "No histogram constraint\n") ;
	fprintf(fp, "Output prefix: %s\n", output_prefix) ;
	fprintf(fp, "-------------------------\n") ;
	fprintf(fp, "iter    time    error\n") ;
	fprintf(fp, "-------------------------\n") ;
	fclose(fp) ;
	
	return 0 ;
}	

int allocate_memory(int flag) {
	algorithm_iterate = malloc(vol * sizeof(float)) ;
	exp_mag = malloc(vol * sizeof(float)) ;
	
	if (flag == 1) {
		obs_mag = malloc(vol * sizeof(float)) ;
		bragg_calc = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
		algorithm_p1 = malloc(vol * sizeof(float)) ;
		algorithm_p2 = malloc(vol * sizeof(float)) ;
		algorithm_r1 = malloc(vol * sizeof(float)) ;
		if (algorithm_beta != 1.)
			algorithm_r2 = malloc(vol * sizeof(float)) ; // for beta != 1
		supp_loc = malloc(vol / 8 * sizeof(long)) ;
	}
	
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	support = malloc(vol * sizeof(uint8_t)) ;
	
	return 0 ;
}

int parse_intens(char *fname, float scale) {
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
				= sqrt(intens[i*s*s + j*s + k]) * scale ;
		else
			obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= intens[i*s*s + j*s + k] ;
	
	free(intens) ;
	
	return 0 ;
}

int parse_bragg(char *fname, double braggqmax) {
	long x, y, z, c = size/2 ;
	double distsq, c2 = c*c ;
	
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
		distsq = ((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) / c2 ;
		
		if (distsq < braggqmax*braggqmax) 
			// Move (q=0) from center to corner
			bragg_calc[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] 
			 = bragg_temp[x*size*size + y*size + z] 
			   * cexpf(-I * 2. * M_PI * (x-c+y-c+z-c) * c / size) ;
		else
			bragg_calc[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] = FLT_MAX ;
	}
	
	fftwf_free(bragg_temp) ;
	
	return 0 ;
}

int parse_support(char *fname) {
	long i ;
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found.\n", fname) ;
		return 1 ;
	}
	fread(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	num_supp = 0 ;
	for (i = 0 ; i < vol ; ++i)
	if (support[i])
		supp_loc[num_supp++] = i ;
	
	fprintf(stderr, "num_supp = %ld\n", num_supp) ;
	supp_index = malloc(num_supp * sizeof(long)) ;
	supp_val = malloc(num_supp * sizeof(float)) ;
	
	return 0 ;
}

void create_plans(char *fname) {
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Measuring plans\n") ;
		forward_plan = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse_plan = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
		
		fp = fopen(fname, "wb") ;
		fftwf_export_wisdom_to_file(fp) ;
		fclose(fp) ;
	
		fprintf(stderr, "Created plans\n") ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		forward_plan = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		inverse_plan = fftwf_plan_dft_3d(size, size, size, fdensity, rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}
}

int gen_input(char *fname, int flag) {
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		if (flag == 0) {
			fprintf(stderr, "Random start\n") ;
			init_model(algorithm_iterate) ;
		}
		else {
			fprintf(stderr, "Cannot find input %s\n", fname) ;
			return 1 ;
		}
	}
	else {
		fprintf(stderr, "Starting from %s\n", fname) ;
		fread(algorithm_iterate, sizeof(float), vol, fp) ;
		fclose(fp) ;
	}
	
	return 0 ;
}

int read_histogram(char *fname, long num) {
	long m, i, num_hist ;
	double frac, sum_hist ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Cannot find histogram file %s\n", fname) ;
		return 1 ;
	}
	fscanf(fp, "%ld", &num_hist) ;
	double *hist = malloc(num_hist * sizeof(double)) ;
	float *val = malloc(num_hist * sizeof(float)) ;
	double *cdf = malloc((num_hist+1) * sizeof(double)) ;
	for (i = 0 ; i < num_hist ; ++i)
		fscanf(fp, "%f %lf", &val[i], &hist[i]) ;
	fclose(fp) ;
	
	cdf[0] = 0. ;
	cdf[1] = hist[0] ;
	sum_hist = hist[0] ;
	for (i = 2 ; i <= num_hist ; ++i) {
		sum_hist += hist[i-1] ;
		cdf[i] = hist[i-1] + cdf[i-1] ;
	}
	for (i = 0 ; i <= num_hist ; ++i)
		cdf[i] /= sum_hist ;
	
	inverse_cdf = malloc(num * sizeof(float)) ;
	inverse_cdf[0] = val[0] ;
	
	i = 1 ;
	for (m = 1 ; m < num ; ++m) {
		frac = (double) m / num ;
		if (i == num_hist - 1)
			inverse_cdf[m] = val[i] ;
		else if (frac <= cdf[i])
			inverse_cdf[m] = val[i-1] + (frac-cdf[i-1]) * (val[i] - val[i-1]) / (cdf[i] - cdf[i-1]) ;
		else {
			i++ ;
			while (frac > cdf[i])
				i++ ;
			inverse_cdf[m] = val[i-1] ;
		}
	}
	
	
	fp = fopen("data/inverse_cdf.raw", "wb") ;
	fwrite(inverse_cdf, sizeof(float), num, fp) ;
	fclose(fp) ;
	
	free(hist) ;
	free(val) ;
	free(cdf) ;
	
	return 0 ;
}
