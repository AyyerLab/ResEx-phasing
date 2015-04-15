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
	
	float blur = 1.0f ;
	char fname[999] ;
	long x, y, z, s = size, c = s/2 ;
	
	fprintf(stderr, "num_supp = %ld\n", num_supp) ; 
	apply_shrinkwrap(iterate, blur, 0.2f) ;
	fprintf(stderr, "Applied shrinkwrap (%.3f). num_supp = %ld\n", blur, num_supp) ;
	
	for (iter = 1 ; iter <= num_iter ; ++iter) {
		gettimeofday(&t1, NULL) ;
		
		error = diffmap(iterate) ;
		
		fprintf(stderr, "\rFinished %d/%d iterations. ", iter, num_iter) ;
		
		if (iter % 50 > 30 || iter%50 == 0) {
			fprintf(stderr, "Now averaging. ") ;
			average_model(p1, average) ;
		}
		
		gettimeofday(&t2, NULL) ;
		fp = fopen("PHASING.log", "a") ;
		fprintf(fp, "%.4d\t%.2f s\t%f\n", 
			iter,
			(double)(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.,
			error) ;
		fclose(fp) ;
		
		if (iter % 50 == 0) {
			blur *= 0.95 ;
			for (i = 0 ; i < vol ; ++i)
				average[i] /= 20. ;
			
			apply_shrinkwrap(average, blur, 0.2f) ;
			fprintf(stderr, "Applied shrinkwrap (%.3f). num_supp = %ld\n", blur, num_supp) ;
			
			sprintf(fname, "data/recon_%d.raw", iter) ;
			fp = fopen(fname, "wb") ;
			fwrite(average, sizeof(float), vol, fp) ;
			fclose(fp) ;
			
			for (i = 0 ; i < vol ; ++i)
				rdensity[i] = average[i] ;
			fftwf_execute(forward_cont) ;
			for (x = 0 ; x < s ; ++x)
			for (y = 0 ; y < s ; ++y)
			for (z = 0 ; z < s ; ++z)
				average[((x+c)%s)*s*s + ((y+c)%s)*s + ((z+c)%s)]
				 = pow(cabsf(fdensity[x*s*s + y*s + z]), 2.) ;
			
			sprintf(fname, "data/frecon_%d.raw", iter) ;
			fp = fopen(fname, "wb") ;
			fwrite(average, sizeof(float), vol, fp) ;
			fclose(fp) ;
			
//			if (iter < num_iter - 1)
				memset(average, 0, vol*sizeof(float)) ;
		}
		
/*		if (iter % 10 == 0) {
			fp = fopen("data/interm.raw", "w") ;
			fwrite(p1, sizeof(float), vol, fp) ;
			fclose(fp) ;
		}
*/	}
	
	fprintf(stderr, "\n") ;
//	fprintf(stderr, "\nCalculating prtf and writing to file.\n") ;
	
//	gen_prtf(average) ;
	
	return 0 ;
}
