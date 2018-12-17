#include "../../src/fft.h"
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol ;
	float r, braggqmax, factor ;
	float *model1 ;
	FILE *fp ;
	char fname[999] ;
	struct fft_data fft ;
	
	if (argc < 4) {
		fprintf(stderr, "Boost Cont: Top hat high pass filter electron density model\n") ;
		fprintf(stderr, "---------------------------------------------------\n") ;
		fprintf(stderr, "Raises magnitudes outside bragg_qmax radius by a constant factor\n") ;
		fprintf(stderr, "\nUsage: %s <model1> <bragg_qmax> <factor>\n", argv[0]) ;
		fprintf(stderr, "\nOutput: <model1>-boost.raw\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	braggqmax = atof(argv[2]) ;
	factor = atof(argv[3]) ;
	sprintf(fname, "%s-boost.raw", remove_ext(argv[1])) ;
	vol = size*size*size ;
	c = size / 2 ;
	
	// Initialize fft
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Parse model
	model1 = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model1, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed model1\n") ;
	
	// Calculate Fourier transform
	for (x = 0 ; x < vol ; ++x)
		fft.rdensity[x] = model1[x] ;
	fft_forward(&fft) ;
	
	// Boost continuous transform
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		r = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		
		if (r > braggqmax * c)
			fft.fdensity[x*size*size + y*size + z] *= factor / vol ;
	}
	
	// Inverse Fourier transform
	fft_inverse(&fft) ;
	
	for (x = 0 ; x < vol ; ++x)
		model1[x] = crealf(fft.rdensity[x]) ;
	
	// Write to file
	fprintf(stderr, "Writing boosted model to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model1, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model1) ;
	fft_free(&fft) ;
	
	return 0 ;
}
