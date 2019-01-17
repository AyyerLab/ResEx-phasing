#include "../../src/fft.h"
#include "../../src/utils.h"
#include "../../src/map.h"

int main(int argc, char* argv[]) {
	long size, vol, num_supp ;
	float blur, thresh ;
	int8_t *support ;
	char fname[1024] ;
	struct fft_data fft ;
	struct ccp4_map map = {0} ;
	
	if (argc < 4) {
		fprintf(stderr, "Create Support: Create support mask from density\n") ;
		fprintf(stderr, "------------------------------------------------\n") ;
		fprintf(stderr, "Gaussian blurring by given width and then thresholding\n") ;
		fprintf(stderr, "\nUsage: %s <density_fname> <blur> <threshold>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "\nOutput: <density_fname>.supp (if <out_fname> not given)\n") ;
		return 1 ;
	}
	blur = atof(argv[2]) ;
	thresh = atof(argv[3]) ;
	
	// Parse density
	parse_map(argv[1], &map) ;
	if (map.f32_data == NULL) {
		fprintf(stderr, "Need float data in input density\n") ;
		return 1 ;
	}
	if (map.header.nx != map.header.ny || map.header.nx != map.header.nz) {
		fprintf(stderr, "Need cubic volume of data (input shape: %d, %d, %d)\n", map.header.nx, map.header.ny, map.header.nz) ;
		return 1 ;
	}
	size = map.header.nx ;
	vol = size*size*size ;
	fprintf(stderr, "Parsed density\n") ;
	
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Apply Shrinkwrap to create support
	support = calloc(vol, sizeof(int8_t)) ;
	num_supp = fft_apply_shrinkwrap(&fft, map.f32_data, blur, thresh, support, NULL) ;
	fprintf(stderr, "Calculated support with %ld voxels\n", num_supp) ;
	
	// Write to file
	if (argc > 4) {
		strcpy(fname, argv[4]) ;
	}
	else {
		sprintf(fname, "%s_supp.ccp4", remove_ext(argv[1])) ;
	}
	int sizes[3] = {size, size, size} ;
	float vsizes[3] = {1.f, 1.f, 1.f} ;
	char label[800] ;
	sprintf(label, "ResEx-phasing:create_support %s %.3e %.3e\n", extract_fname(argv[1]), blur, thresh) ;
	int flipped = (map.header.mapc == 3) ;
	fprintf(stderr, "Writing support to %s\n", fname) ;
	save_mask_as_map(fname, support, sizes, vsizes, label, flipped) ;
	
	// Free memory
	free_map(&map) ;
	free(support) ;
	
	return 0 ;
}
