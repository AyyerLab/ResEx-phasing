#include "../../src/fft.h"
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	long i, size, vol ;
	float *temp, r_max = -1 ;
	char fname[1024], label[800] ;
	struct fft_data fft ;
	struct ccp4_map map = {0} ;
	
	if (argc < 2) {
		fprintf(stderr, "Gen Dens: Inverse fourier transform complex amplitudes to get density\n") ;
		fprintf(stderr, "---------------------------------------------------------------------\n") ;
		fprintf(stderr, "Can supply cutoff radius in voxels to truncate resolution of output\n") ;
		fprintf(stderr, "\nUsage: %s <fmodel_name>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "Second option: <r_max> cutoff radius\n") ;
		fprintf(stderr, "\nOutput: <cpx_fname>-dens.ccp4 (if <out_fname> not given)\n") ;
		return 1 ;
	}
	if (argc > 2)
		strcpy(fname, argv[2]) ;
	else
		sprintf(fname, "%s-dens.raw", remove_ext(argv[1])) ;
	
	if (argc > 3) {
		r_max = atof(argv[3]) ;
		fprintf(stderr, "Truncated to radius = %f\n", r_max) ;
	}
	
	// Read complex Fourier amplitudes
	if (parse_map(argv[1], &map))
		return 1 ;
	if (map.c64_data == NULL) {
		fprintf(stderr, "Need float complex data in input fmodel\n") ;
		return 1 ;
	}
	if (map.header.nx != map.header.ny || map.header.nx != map.header.nz) {
		fprintf(stderr, "Need cubic volume of data (input shape: %d, %d, %d)\n", map.header.nx, map.header.ny, map.header.nz) ;
		return 1 ;
	}
	size = map.header.nx ;
	vol = size*size*size ;
	fprintf(stderr, "Parsed fdensity ( size = %ld )\n", size) ;
	
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Shift coordinates such that origin is at the corner
	// Also truncate to radius if specified
	fft_shift_complex(&fft, fft.fdensity, map.c64_data, r_max) ;
	
	// Do inverse Fourier transform
	fft_inverse(&fft) ;
	
	// Shift real space density such that it is in the center of the cube
	temp = calloc(vol, sizeof(float)) ;
	fft_ishift_real(&fft, temp, fft.rdensity, -1.) ;
	for (i = 0 ; i < vol ; ++i)
		temp[i] /= vol ;
	
	fprintf(stderr, "Writing density to %s\n", fname) ;
	int sizes[3] = {size, size, size} ;
	float vsizes[3] = {1.f, 1.f, 1.f} ;
	sprintf(label, "ResEx-phasing:gen_dens %s\n", extract_fname(argv[1])) ;
	int flipped = (map.header.mapc == 3) ;
	save_vol_as_map(fname, temp, sizes, vsizes, label, flipped) ;
	
	free_map(&map) ;
	free(temp) ;
	fft_free(&fft) ;
	
	return 0 ;
}

