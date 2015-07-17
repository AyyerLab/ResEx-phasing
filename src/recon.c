#include "brcont.h"

int main(int argc, char *argv[]) {
	int i, num_iter, start_ave ;
	double error ;
	float *average ;
	struct timeval t1, t2 ;
	FILE *fp ;
	
	if (argc < 3) {
//		fprintf(stderr, "Format: %s <num_iter> <start_ave>\n", argv[0]) ;
		fprintf(stderr, "Format: %s <num_iter> <start_ER>\n", argv[0]) ;
		return 1 ;
	}
	num_iter = atoi(argv[1]) ;
	start_ave = atoi(argv[2]) ;
	
	if (start_ave > num_iter) {
		fprintf(stderr, "start_ave > num_iter, printing last iteration\n") ;
		start_ave = num_iter ;
	}
	
	if (setup())
		return 2 ;
	
	omp_set_num_threads(32) ;
	
	average = calloc(vol, sizeof(float)) ;
	
	fp = fopen("PHASING.log", "w") ;
	fprintf(fp, "iter\ttime  \terror\n") ;
	fprintf(fp, "-------------------------\n") ;
	fclose(fp) ;
	
	fprintf(stderr, "num_supp = %ld\n", num_supp) ; 
	
	for (iter = 1 ; iter <= num_iter ; ++iter) {
		gettimeofday(&t1, NULL) ;
		
//		if (iter < start_ave)
			error = diffmap(iterate) ;
//		else {
//			fprintf(stderr, "Doing error-reduction. ") ;
//			error = error_red(iterate) ;
//		}
		
		fprintf(stderr, "\rFinished %d/%d iterations. ", iter, num_iter) ;
		
		if (iter >= start_ave) {
			fprintf(stderr, "Now averaging. ") ;
			average_model(p1, average) ;
		}

		//if (iter == 100 || iter == 150)
		//	apply_shrinkwrap(p1, 1.5, 0.2) ;
		
		gettimeofday(&t2, NULL) ;
		fp = fopen("PHASING.log", "a") ;
		fprintf(fp, "%.4d\t%.2f s\t%f\n", 
			iter,
			(double)(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.,
			error) ;
		fclose(fp) ;
		
		if (iter % 20 == 0) {
			fp = fopen("data/interm.raw", "w") ;
			fwrite(p1, sizeof(float), vol, fp) ;
			fclose(fp) ;
		}
	}
	
	fprintf(stderr, "\nCalculating prtf and writing to file.\n") ;
	
	for (i = 0 ; i < vol ; ++i)
		average[i] /= (num_iter - start_ave + 1) ;
//		average[i] = iterate[i] ;
	
	fp = fopen("data/recon.raw", "wb") ;
	fwrite(average, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	gen_prtf(average) ;
	
	return 0 ;
}
