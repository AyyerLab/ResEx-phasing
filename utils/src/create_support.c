#include "../../src/fft.h"
#include "../../src/utils.h"

int main(int argc, char* argv[]) {
	long size, vol, num_supp ;
	float blur, thresh, *density ;
	uint8_t *support ;
	char fname[1024] ;
	FILE *fp ;
	struct fft_data fft ;
	
	if (argc < 4) {
		fprintf(stderr, "Create Support: Create support mask from density\n") ;
		fprintf(stderr, "------------------------------------------------\n") ;
		fprintf(stderr, "Gaussian blurring by given width and then thresholding\n") ;
		fprintf(stderr, "\nUsage: %s <density_fname> <blur> <threshold>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname>\n") ;
		fprintf(stderr, "\nOutput: <density_fname>.supp (if <out_fname> not given)\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	blur = atof(argv[2]) ;
	thresh = atof(argv[3]) ;
	
	vol = size*size*size ;
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Parse density
	density = malloc(vol * sizeof(float)) ;
	support = calloc(vol, sizeof(uint8_t)) ;
	fp = fopen(argv[1], "rb") ;
	fread(density, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed density\n") ;
	
	// Apply Shrinkwrap to create support
	num_supp = fft_apply_shrinkwrap(&fft, density, blur, thresh, support, NULL) ;
	fprintf(stderr, "Calculated support with %ld voxels\n", num_supp) ;
	
	// Write to file
	if (argc > 4) {
		strcpy(fname, argv[4]) ;
	}
	else {
		sprintf(fname, "%s.supp", remove_ext(argv[1])) ;
	}
	fprintf(stderr, "Writing support to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(density) ;
	free(support) ;
	
	return 0 ;
}
