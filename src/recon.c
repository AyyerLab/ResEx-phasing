#include "brcont.h"

int main(int argc, char *argv[]) {
	int c, i, num_iter = -1, start_ave = -1 ;
	double error ;
	float *average ;
	struct timeval t1, t2 ;
	FILE *fp ;
	char fname[999], config_fname[999] ;
	
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
		else {
			num_iter = atoi(argv[optind++]) ;
			if (optind >= argc) {
				fprintf(stderr, "Missing start_ave\n") ;
				break ;
			}
			start_ave = atoi(argv[optind++]) ;
		}
	}
	
	if (num_iter == -1)
		fprintf(stderr, "Missing num_iter\n") ;
	
	if (num_iter == -1 || start_ave == -1) {
		fprintf(stderr, "Format: %s [-c config_fname] num_iter start_ave\n", argv[0]) ;
		fprintf(stderr, "Default: -c config.ini\n") ;
		return 1 ;
	}
	
	if (start_ave > num_iter) {
		fprintf(stderr, "start_ave > num_iter, printing last iteration\n") ;
		start_ave = num_iter ;
	}
	
	if (setup(config_fname))
		return 2 ;
	
	omp_set_num_threads(omp_get_max_threads()) ;
	
	average = calloc(vol, sizeof(float)) ;
	
	sprintf(fname, "%s-log.dat", output_prefix) ;
	fp = fopen(fname, "w") ;
	fprintf(fp, "iter\ttime  \terror\n") ;
	fprintf(fp, "-------------------------\n") ;
	fclose(fp) ;
	
	fprintf(stderr, "num_supp = %ld\n", num_supp) ; 
	
	for (iter = 1 ; iter <= num_iter ; ++iter) {
		gettimeofday(&t1, NULL) ;
		
//		if (iter < start_ave)
			error = diffmap(iterate) ;
//			error = modified_hio(iterate) ;
//		else {
//			fprintf(stderr, "Doing error-reduction. ") ;
//			error = error_red(iterate) ;
//		}
		
		fprintf(stderr, "\rFinished %d/%d iterations. ", iter, num_iter) ;
		
		if (iter >= start_ave) {
//			fprintf(stderr, "Doing ER. ") ;
			fprintf(stderr, "Now averaging. ") ;
			average_model(p1, average) ;
		}
		
		gettimeofday(&t2, NULL) ;
		sprintf(fname, "%s-log.dat", output_prefix) ;
		fp = fopen(fname, "a") ;
		fprintf(fp, "%.4d\t%.2f s\t%f\n", 
			iter,
			(double)(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.,
			error) ;
		fclose(fp) ;
	}
	
	fprintf(stderr, "\nCalculating prtf and writing to file.\n") ;
	
	for (i = 0 ; i < vol ; ++i)
		average[i] /= (num_iter - start_ave + 1) ;
//		average[i] = iterate[i] ;
	
	sprintf(fname, "%s-recon.raw", output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(average, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	gen_prtf(average) ;
	
	return 0 ;
}
