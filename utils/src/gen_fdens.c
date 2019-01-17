#include "../../src/fft.h"
#include "../../src/utils.h"
#include "../../src/volume.h"
#include "../../src/map.h"

int main(int argc, char *argv[]) {
	long x, size, vol ;
	float *temp ;
	char fname[1024], symfname[1024], label[800] ;
	struct fft_data fft ;
	struct volume_data volume ;
	struct ccp4_map map = {0} ;
	strcpy(volume.point_group, "222") ;
	
	if (argc < 2) {
		fprintf(stderr, "Gen FDens: Fourier transform electron densities\n") ;
		fprintf(stderr, "-----------------------------------------------\n") ;
		fprintf(stderr, "Densities assumed to cube of 32-bit floats\n") ;
		fprintf(stderr, "Also produces output symmetrized by point group (default '1')\n") ;
		fprintf(stderr, "\nUsage: %s <model_fname>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname> <point_group>\n") ;
		fprintf(stderr, "\nOutput: <model_fname>-fdens.ccp4 (if <out_fname> not given)\n") ;
		fprintf(stderr, "\t<out_fname>-sym.ccp4\n") ;
		return 1 ;
	}
	if (argc > 2) {
		strcpy(fname, argv[2]) ;
		sprintf(symfname, "%s-sym.ccp4", remove_ext(fname)) ;
	}
	else {
		sprintf(fname, "%s-fdens.ccp4", remove_ext(argv[1])) ;
		sprintf(symfname, "%s-fdens-sym.ccp4", remove_ext(argv[1])) ;
	}
	if (argc > 3)
		strcpy(volume.point_group, argv[3]) ;
	
	// Parse real-space density
	parse_map(argv[1], &map) ;
	if (map.f32_data == NULL) {
		fprintf(stderr, "Need float data in input model\n") ;
		return 1 ;
	}
	if (map.header.nx != map.header.ny || map.header.nx != map.header.nz) {
		fprintf(stderr, "Need cubic volume of data (input shape: %d, %d, %d)\n", map.header.nx, map.header.ny, map.header.nz) ;
		return 1 ;
	}
	size = map.header.nx ;
	vol = size*size*size ;
	fprintf(stderr, "Parsed density ( size = %ld )\n", size) ;
	
	volume_init(&volume, size) ;
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Translate array such that molecule is around the origin
	for (x = 0 ; x < vol ; ++x)
		fft.fdensity[x] = map.f32_data[x] ;
	fft_shift_complex(&fft, fft.rdensity, fft.fdensity, -1.) ;
	
	// Apply Fourier transform
	fft_forward(&fft) ;
	fprintf(stderr, "Fourier transformed density\n") ;
	
	// Translate array to put the origin at (c, c, c)
	fft_ishift_complex(&fft, fft.rdensity, fft.fdensity, -1.) ;
	
	// Write complex array to file
	fprintf(stderr, "Writing fdensity to %s\n", fname) ;
	int sizes[3] = {size, size, size} ;
	float vsizes[3] = {-1.f, -1.f, -1.f} ; // Negative to indicate Fourier space
	sprintf(label, "ResEx-phasing:gen_fdens %s\n", extract_fname(argv[1])) ;
	int flipped = (map.header.mapc == 3) ;
	save_cpx_as_map(fname, fft.rdensity, sizes, vsizes, label, flipped) ;
	
	// Symmetrize intensities
	temp = calloc(vol, sizeof(float)) ;
	volume_symmetrize_centered(&volume, fft.rdensity, temp) ; 
	
	// Write symmetrized intensities to file
	fprintf(stderr, "Writing symmetrized intensity to %s\n", symfname) ;
	sprintf(label, "ResEx-phasing:gen_fdens(sym) %s %s\n", extract_fname(argv[1]), volume.point_group) ;
	save_vol_as_map(symfname, temp, sizes, vsizes, label, flipped) ;
	
	// Free memory
	free_map(&map) ;
	free(temp) ;
	fft_free(&fft) ;
	
	return 0 ;
}

