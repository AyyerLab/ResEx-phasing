#include "../../src/fft.h"
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long i, size, vol ;
	float *temp, r_max = -1 ;
	FILE *fp ;
	char fname[999] ;
	struct fft_data fft ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <cpx_fmodel>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "Second option: <r_max> cutoff radius\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float complex)) ;
	vol = size*size*size ;
	
	if (argc > 2)
		strcpy(fname, argv[2]) ;
	else
		sprintf(fname, "%s-dens.raw", remove_ext(argv[1])) ;
	
	if (argc > 3) {
		r_max = atof(argv[3]) ;
		fprintf(stderr, "Truncated to radius = %f\n", r_max) ;
	}
	
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	temp = malloc(vol * sizeof(float)) ;
	
	// Read complex Fourier amplitudes
	fp = fopen(argv[1], "rb") ;
	fread(fft.rdensity, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Shift coordinates such that origin is at the corner
	// Also truncate to radius if specified
	fft_shift_complex(&fft, fft.fdensity, fft.rdensity, r_max) ;
	
	// Do inverse Fourier transform
	fft_inverse(&fft) ;
	
	// Shift real space density such that it is in the center of the cube
	fft_ishift_real(&fft, temp, fft.rdensity, -1.) ;
	for (i = 0 ; i < vol ; ++i)
		temp[i] /= vol ;
	
	fprintf(stderr, "Writing density to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(temp) ;
	fft_free(&fft) ;
	
	return 0 ;
}

