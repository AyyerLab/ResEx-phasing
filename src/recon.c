#include "brcont.h"

int parse_args(int, char*[], char*) ;

int main(int argc, char *argv[]) {
	long i, num_vox ;
	double error ;
	float *average_p1, *average_p2 ;
	struct timeval t1, t2 ;
	FILE *fp ;
	char algo[8], fname[999], config_fname[999] ;
	
	if (parse_args(argc, argv, config_fname))
		return 1 ;
	
	if (setup(config_fname))
		return 2 ;
	
	omp_set_num_threads(omp_get_max_threads()) ;
	if (do_bg_fitting)
		num_vox = 2 * vol ;
	else
		num_vox = vol ;
	
	average_p1 = calloc(num_vox, sizeof(float)) ;
	average_p2 = calloc(num_vox, sizeof(float)) ;
	sprintf(fname, "%s-slices", output_prefix) ;
	mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	sprintf(fname, "%s-fslices", output_prefix) ;
	mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	if (do_bg_fitting) {
		sprintf(fname, "%s-radavg", output_prefix) ;
		mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	}
	if (do_local_variation) {
		sprintf(fname, "%s-support", output_prefix) ;
		mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	}
	
	for (iter = 1 ; iter <= num_iter+num_avg_iter ; ++iter) {
		gettimeofday(&t1, NULL) ;
		
		if (iter <= num_iter)
			strcpy(algo, algorithms[iter-1]) ;
		else
			strcpy(algo, avg_algorithms[iter-num_iter-1]) ;
		
		if (strcmp(algo, "DM") == 0)
			error = DM_algorithm(algorithm_iterate) ;
		else if (strcmp(algo, "HIO") == 0)
			error = HIO_algorithm(algorithm_iterate) ;
		else if (strcmp(algo, "RAAR") == 0)
			error = RAAR_algorithm(algorithm_iterate) ;
		else if (strcmp(algo, "mod-DM") == 0)
			error = mod_DM_algorithm(algorithm_iterate) ;
		else if (strcmp(algo, "ER") == 0)
			error = ER_algorithm(algorithm_iterate) ;
		else {
			fprintf(stderr, "Could not understand algorithm name: %s\n", algo) ;
			return 1 ;
		}
		
		if (iter > num_iter) {
			average_model(algorithm_p1, average_p1) ;
			average_model(algorithm_p2, average_p2) ;
		}
		
		gettimeofday(&t2, NULL) ;
		sprintf(fname, "%s-log.dat", output_prefix) ;
		fp = fopen(fname, "a") ;
		fprintf(fp,
		        "%.4d\t%.2f s\t%f\n",
		        iter,
			    (double)(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.,
			    error) ;
		fclose(fp) ;
		
		if (iter == num_iter) {
			sprintf(fname, "%s-last.raw", output_prefix) ;
			fp = fopen(fname, "w") ;
			fwrite(algorithm_iterate, sizeof(float), num_vox, fp) ;
			fclose(fp) ;
		}
		
		if (iter%1 == 0) {
			sprintf(fname, "%s-slices/%.4d.raw", output_prefix, iter) ;
			dump_slices(algorithm_p1, fname, 0) ;
			sprintf(fname, "%s-fslices/%.4d.raw", output_prefix, iter) ;
			dump_slices(exp_mag, fname, 1) ;
			if (do_local_variation) {
				sprintf(fname, "%s-support/%.4d.supp", output_prefix, iter) ;
				dump_support_slices(support, fname) ;
			}
			if (do_bg_fitting) {
				sprintf(fname, "%s-radavg/%.4d.raw", output_prefix, iter) ;
				fp = fopen(fname, "wb") ;
				fwrite(radavg, sizeof(double), size/2, fp) ;
				fclose(fp) ;
			}
		}

//		if (iter > 250 && iter % 20 == 0)
//			num_supp = (0.99 * num_supp) ;
		
		fprintf(stderr, "\rFinished %d/%d iterations. ", iter, num_iter+num_avg_iter) ;
//		fprintf(stderr, "\rFinished %d/%d iterations. num_supp = %ld. ", iter, num_iter+num_avg_iter, num_supp) ;
		if (iter > num_iter)
			fprintf(stderr, "Now averaging") ;
	}
	
	fprintf(stderr, "\nCalculating prtf and writing to file.\n") ;
	
	if (num_avg_iter > 0) {
		for (i = 0 ; i < num_vox ; ++i) {
			average_p1[i] /= num_avg_iter ;
			average_p2[i] /= num_avg_iter ;
		}
	}
	else {
		for (i = 0 ; i < num_vox ; ++i) {
			average_p1[i] = algorithm_p1[i] ;
			average_p2[i] = algorithm_p2[i] ;
		}
	}
	
	sprintf(fname, "%s-pf.raw", output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(average_p1, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	sprintf(fname, "%s-pd.raw", output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(average_p2, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	if (do_bg_fitting) {
		sprintf(fname, "%s-bg.raw", output_prefix) ;
		fp = fopen(fname, "wb") ;
		fwrite(&(algorithm_p1[vol]), sizeof(float), vol, fp) ;
		fclose(fp) ;
		
		sprintf(fname, "%s-radavg.raw", output_prefix) ;
		fp = fopen(fname, "wb") ;
		fwrite(radavg, sizeof(float), size/2, fp) ;
		fclose(fp) ;
	}
	
	if (do_local_variation) {
		sprintf(fname, "%s-supp.supp", output_prefix) ;
		fp = fopen(fname, "wb") ;
		fwrite(support, sizeof(uint8_t), vol, fp) ;
		fclose(fp) ;
	}
	
	gen_prtf(average_p2) ;
	
	free(average_p1) ;
	free(average_p2) ;
	
	return 0 ;
}

int parse_args(int argc, char *argv[], char *config_fname) {
	int c ;
	extern char *optarg ;
	extern int optind ;
	
	strcpy(config_fname, "config.ini") ;
	
	while (optind < argc) {
		if ((c = getopt(argc, argv, "c:")) != -1) {
			switch (c) {
				case 'c':
					strcpy(config_fname, optarg) ;
					break ;
			}
		}
	}
	
	return 0 ;
}

