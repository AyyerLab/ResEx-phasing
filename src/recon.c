#include "algorithm.h"

int parse_args(int argc, char *argv[], char *config_fname) {
	extern char *optarg ;
	extern int optind ;
	int testing_mode = 0 ;
	
	strcpy(config_fname, "config.ini") ;
	
	while (optind < argc) {
		int c ;
		if ((c = getopt(argc, argv, "c:Th")) != -1) {
			switch (c) {
				case 'c':
					strcpy(config_fname, optarg) ;
					break ;
				case 'T':
					testing_mode = 1 ;
					break ;
				case 'h':
					fprintf(stderr, "Format: %s [-c <config.ini>] [-T] [-h]\n", argv[0]) ;
					return -1 ;
					break ;
				case '?':
					fprintf(stderr, "Invalid command line option\n") ;
					return -1 ;
			}
		}
	}
	
	return testing_mode ;
}

int setup(struct algorithm_data *self, char *config_fname, int fixed_seed) {
	char line[1024] ;
	char input_fname[1024], hist_fname[1024] ;
	char intens_fname[1024], bragg_fname[1024] ;
	char support_fname[1024] ;
	char inputbg_fname[1024], quat_fname[1024] ;
	char algorithm_string[8192], avg_algorithm_string[8192] ;
	float bragg_qmax = 0., sigma = 0. ;
	float scale_factor = 0. ;
	int num_threads = -1, num_div = -1 ;
	long size = 0 ;
	self->volume = malloc(sizeof(struct volume_data)) ;
	self->input = malloc(sizeof(struct input_data)) ;
	self->fft = malloc(sizeof(struct fft_data)) ;
	self->quat = malloc(sizeof(struct rotation)) ;
	struct volume_data *volume = self->volume ;
	struct input_data *input = self->input ;
	struct fft_data *fft = self->fft ;
	struct rotation *quat = self->quat ;
	
	self->output_prefix[0] = '\0' ;
	strcpy(volume->point_group, "222") ;
	self->do_blurring = 0 ;
	self->do_histogram = 0 ;
	self->do_positivity = 0 ;
	self->do_local_variation = 0 ;
	self->do_bg_fitting = 0 ;
	self->do_normalize_prtf = 0 ;
	strcpy(algorithm_string, "DM") ;
	strcpy(avg_algorithm_string, "Avg") ;
	
	FILE *fp = fopen(config_fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Config file, %s, not found.\n", config_fname) ;
		return 1 ;
	}
	while (fgets(line, 1024, fp) != NULL) {
		char *token = strtok(line, " =") ;
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
			strcpy(volume->point_group, strtok(NULL, " =\n")) ;
		// Files
		else if (strcmp(token, "intens_fname") == 0)
			strcpy(intens_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "bragg_fname") == 0)
			strcpy(bragg_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "input_fname") == 0)
			strcpy(input_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "inputbg_fname") == 0)
			strcpy(inputbg_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "support_fname") == 0)
			strcpy(support_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "output_prefix") == 0)
			strcpy(self->output_prefix, strtok(NULL, " =\n")) ;
		// Algorithm
		else if (strcmp(token, "algorithm") == 0)
			strcpy(algorithm_string, strtok(NULL, "=\n")) ;
		else if (strcmp(token, "avg_algorithm") == 0)
			strcpy(avg_algorithm_string, strtok(NULL, "=\n")) ;
		else if (strcmp(token, "beta") == 0)
			self->beta = atof(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "bg_fitting") == 0)
			self->do_bg_fitting = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "blurring") == 0)
			self->do_blurring = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "histogram") == 0)
			self->do_histogram = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "local_variation") == 0)
			self->do_local_variation = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "positivity") == 0)
			self->do_positivity = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "normalize_prtf") == 0)
			self->do_normalize_prtf = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "quat_fname") == 0)
			strcpy(quat_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "num_div") == 0)
			num_div = atoi(strtok(NULL, " =\n")) ;
		else if (strcmp(token, "hist_fname") == 0)
			strcpy(hist_fname, strtok(NULL, " =\n")) ;
		else if (strcmp(token, "sigma_deg") == 0)
			sigma = atof(strtok(NULL, " =\n")) ;
		else
			fprintf(stderr, "Unable to recognise option: \"%s\"\n", token) ;
	}
	fclose(fp) ;
	
	if (size == 0) {
		fprintf(stderr, "Need nonzero size in config file\n") ;
		return 1 ;
	}
	else if (size%2 == 0) {
		fprintf(stderr, "size = %ld is even. An odd number is preferred.\n", size) ;
	}
	if (strcmp(volume->point_group, "222")*strcmp(volume->point_group, "4")*strcmp(volume->point_group, "1") != 0) {
		fprintf(stderr, "Point group needs to be either 1, 222 or 4 for now\n") ;
		return 1 ;
	}
	if (self->output_prefix[0] == '\0') {
		fprintf(stderr, "Using default output prefix data/output\n") ;
		strcpy(self->output_prefix, "data/output") ;
	}
	
	self->size = size ;
	self->vol = size*size*size ;
	self->num_vox = self->do_bg_fitting ? self->vol * 2 : self->vol ;
	num_threads = (num_threads == -1 ? omp_get_max_threads() : num_threads) ;
	omp_set_num_threads(num_threads) ;
	
	fft_init(fft, size, num_threads) ;
	fft_create_plans(fft) ;
	input_init(input, size) ;
	volume_init(volume, size) ;
	
	if (parse_algorithm_strings(self, algorithm_string, avg_algorithm_string))
		return 1 ;
	algorithm_allocate_memory(self) ;
	if (input_parse_intens(input, intens_fname, scale_factor, self->do_bg_fitting))
		return 1 ;
	if (input_parse_bragg(input, bragg_fname, bragg_qmax))
		return 1 ;
	if (input_parse_support(input, support_fname))
		return 1 ;
	if (self->do_histogram) {
		if (input_read_histogram(input, hist_fname))
			return 1 ;
	}
	if (self->do_blurring) {
		if (quat_fname[0] != '\0' && num_div > 0) {
			fprintf(stderr, "quat_fname and num_div both specified. Pick one.\n") ;
			return 1 ;
		}
		else if (quat_fname[0] != '\0') {
			if (quat_parse(quat, quat_fname))
				return 1 ;
		}
		else if (num_div > 0) {
			if (sigma == 0.) {
				fprintf(stderr, "num_div option also requires a sigma_deg value to subset the quaternion file\n") ;
				return 1 ;
			}
			quat_gen(quat, num_div, sigma) ;
		}
		else {
			fprintf(stderr, "Need either num_div or quat_fname if do_blurring is active\n") ;
			return 1 ;
		}
	}
	
	input_init_iterate(input, input_fname, inputbg_fname, self->iterate, self->do_bg_fitting, fixed_seed) ;
	fft_create_plans(fft) ;
	if (self->do_bg_fitting)
		volume_init_radavg(volume) ;
	
	sprintf(line, "%s-log.dat", self->output_prefix) ;
	fp = fopen(line, "w") ;
	fprintf(fp, "Resolution extension iterative phasing\n") ;
	fprintf(fp, "Data: %s %s\n", bragg_fname, intens_fname) ;
	fprintf(fp, "Support: %s (%ld)\n", support_fname, input->num_supp) ;
	fprintf(fp, "Algorithm: %s with beta = %.2f\n", algorithm_string, self->beta) ;
	fprintf(fp, "Averaging algorithm: %s\n", avg_algorithm_string) ;
	if (self->do_positivity)
		fprintf(fp, "Assuming electron density is positive\n") ;
	if (self->do_histogram)
		fprintf(fp, "Applying histogram constraint: %s\n", hist_fname) ;
	if (self->do_local_variation)
		fprintf(fp, "Updating support using local variation\n") ;
	if (self->do_bg_fitting)
		fprintf(fp, "Fitting spherically symmetric background\n") ;
	if (self->do_blurring)
		fprintf(fp, "Rotationally blurring model with %d orientations\n", quat->num_rot) ;
	if (self->do_normalize_prtf)
		fprintf(fp, "Normalizing output by PRTF\n") ;
	fprintf(fp, "Output prefix: %s\n", self->output_prefix) ;
	fprintf(fp, "-------------------------\n") ;
	fprintf(fp, "iter    time    error\n") ;
	fprintf(fp, "-------------------------\n") ;
	fclose(fp) ;
	
	make_recon_folders(self) ;
	
	return 0 ;
}

int main(int argc, char *argv[]) {
	long i ;
	int iter, fixed_seed ;
	struct timeval t1, t2 ;
	char config_fname[1024] ;
	
	struct algorithm_data algo ;
	
	if ((fixed_seed = parse_args(argc, argv, config_fname)) < 0)
		return 1 ;
	
	if (setup(&algo, config_fname, fixed_seed))
		return 2 ;
	
	for (iter = 1 ; iter <= algo.num_iter+algo.num_avg_iter ; ++iter) {
		gettimeofday(&t1, NULL) ;
		
		float error = run_iteration(&algo, iter) ;
		if (error < 0)
			return 1 ;
		
		if (iter > algo.num_iter) {
			volume_accumulate(algo.p1, algo.average_p1, algo.num_vox) ;
			volume_accumulate(algo.p2, algo.average_p2, algo.num_vox) ;
		}
		
		gettimeofday(&t2, NULL) ;
		save_current(&algo, iter, t1, t2, error) ;
		
		fprintf(stderr, "\rFinished %d/%d iterations. ", iter, algo.num_iter+algo.num_avg_iter) ;
		if (iter > algo.num_iter)
			fprintf(stderr, "Now averaging") ;
	}
	fprintf(stderr, "\nCalculating prtf and writing to file.\n") ;
	
	if (algo.num_avg_iter > 0) {
		for (i = 0 ; i < algo.num_vox ; ++i) {
			algo.average_p1[i] /= algo.num_avg_iter ;
			algo.average_p2[i] /= algo.num_avg_iter ;
		}
	}
	else {
		for (i = 0 ; i < algo.num_vox ; ++i) {
			algo.average_p1[i] = algo.p1[i] ;
			algo.average_p2[i] = algo.p2[i] ;
		}
	}
	
	calc_prtf(&algo, 100) ;
	save_output(&algo) ;
	
	return 0 ;
}

