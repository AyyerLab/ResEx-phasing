#include "../../src/fft.h"
#include "../../src/utils.h"
#include "../../src/volume.h"

int main(int argc, char *argv[]) {
	long x, size, vol ;
	float *temp ;
	FILE *fp ;
	char fname[999], symfname[999] ;
	struct fft_data fft ;
	struct volume_data volume ;
	strcpy(volume.point_group, "222") ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <raw_model>\n", argv[0]) ;
		fprintf(stderr, "Optional: <out_fname> <point_group>\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	if (argc > 2) {
		strcpy(fname, argv[2]) ;
		sprintf(symfname, "%s-sym.raw", remove_ext(fname)) ;
	}
	else {
		sprintf(fname, "%s-fdens.cpx", remove_ext(argv[1])) ;
		sprintf(symfname, "%s-fdens-sym.raw", remove_ext(argv[1])) ;
	}
	if (argc > 3)
		strcpy(volume.point_group, argv[3]) ;
	vol = size*size*size ;
	
	volume_init(&volume, size) ;
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	// Allocate memory
	temp = malloc(vol * sizeof(float)) ;
	
	// Parse real-space density
	fp = fopen(argv[1], "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed density\n") ;
	
	// Translate array such that molecule is around the origin
	for (x = 0 ; x < vol ; ++x)
		fft.fdensity[x] = temp[x] ;
	fft_shift_complex(&fft, fft.rdensity, fft.fdensity, -1.) ;
	
	// Apply Fourier transform
	fft_forward(&fft) ;
	fprintf(stderr, "Fourier transformed density\n") ;
	
	// Translate array to put the origin at (c, c, c)
	fft_ishift_complex(&fft, fft.rdensity, fft.fdensity, -1.) ;
	
	// Write complex array to file
	fprintf(stderr, "Writing fdensity to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(fft.rdensity, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	// Symmetrize intensities
	memset(temp, 0, vol * sizeof(float)) ;
	volume_symmetrize_centered(&volume, fft.rdensity, temp) ; 
	
	// Write symmetrized intensities to file
	fprintf(stderr, "Writing symmetrized intensity to %s\n", symfname) ;
	fp = fopen(symfname, "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Free memory
	free(temp) ;
	fft_free(&fft) ;
	
	return 0 ;
}

