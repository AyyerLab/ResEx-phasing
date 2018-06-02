#include "../../src/utils.h"
#include "../../src/fft.h"

int main(int argc, char *argv[]) {
	long x, y, z, i, size, c, vol ;
	float *recon, *temp ;
	FILE *fp ;
	char fname[999] ;
	struct fft_data fft ;
	
	if (argc < 3) {
		fprintf(stderr, "Averages a set of reconstructions in the results folder\n") ;
		fprintf(stderr, "Format: %s <num1> <num2> ...\n", argv[0]) ;
		return 1 ;
	}
	size = 501 ;
	c = size / 2 ;
	vol = size*size*size ;
	
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	recon = malloc(vol * sizeof(float)) ;
	temp = malloc(vol * sizeof(float)) ;
	
	// Calculate mean real-space density
	for (x = 0 ; x < vol ; ++x)
		fft.rdensity[x] = 0.f ;
	
	for (i = 1 ; i < argc ; ++i) {
		sprintf(fname, "results/output_%.2d.raw", atoi(argv[i])) ;
		fp = fopen(fname, "rb") ;
		fread(recon, sizeof(float), vol, fp) ;
		fclose(fp) ;
		
		for (x = 0 ; x < vol ; ++x)
			fft.rdensity[x] += recon[x] ;
	}
	
	for (x = 0 ; x < vol ; ++x) {
		fft.rdensity[x] /= (argc - 1) ;
		temp[x] = crealf(fft.rdensity[x]) ;
	}
	
	// Write mean real-space density to file
	fp = fopen("results/comb_output.raw", "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate mean Fourier intensity
	fft_forward(&fft) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
		  = powf(cabsf(fft.fdensity[x*size*size + y*size + z]), 2.f) ;
	
	fp = fopen("results/comb_foutput.raw", "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	free(recon) ;
	fft_free(&fft) ;
	
	return 0 ;
}

