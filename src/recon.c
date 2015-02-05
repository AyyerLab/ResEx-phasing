#include "brcont.h"

int main(int argc, char *argv[]) {
	int i, num_iter, start_ave ;
	double error ;
	float *average ;
	struct timeval t1, t2 ;
	FILE *fp ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <num_iter> <start_ave>\n", argv[0]) ;
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
	
	for (iter = 1 ; iter <= num_iter ; ++iter) {
		gettimeofday(&t1, NULL) ;
		
		error = diffmap(iterate) ;
		
		if (iter >= start_ave) {
			average_model(p1, average) ;
		}
		
		gettimeofday(&t2, NULL) ;
		fp = fopen("PHASING.log", "a") ;
		fprintf(fp, "%.4d\t%.2f s\t%f\n", 
			iter,
			(double)(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.,
			error) ;
		fclose(fp) ;
	}
	
	for (i = 0 ; i < vol ; ++i)
		average[i] /= (num_iter - start_ave + 1) ;
	
	gen_prtf(average) ;
	
	fp = fopen("data/recon.raw", "wb") ;
	fwrite(average, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	return 0 ;
}
