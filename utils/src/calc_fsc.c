#include "../../src/utils.h"
#include "../../src/fft.h"

void process_model(char *fname, float complex *model, struct fft_data *fft, float cm[3]) {
	long x, y, z ;
	long size = fft->size, c = size/2, vol = size*size*size ;
	FILE *fp ;
	float *temp = malloc(vol * sizeof(float)) ;
	float cm_denr = 0. ;
	
	// Parse first model
	fp = fopen(fname, "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate center of mass
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		cm[0] += x*temp[x*size*size + y*size + z] ;
		cm[1] += y*temp[x*size*size + y*size + z] ;
		cm[2] += z*temp[x*size*size + y*size + z] ;
		cm_denr += temp[x*size*size + y*size + z] ;
	}
	cm[0] /= cm_denr ; cm[1] /= cm_denr ; cm[2] /= cm_denr ;
	
	// Calculate Fourier transform
	for (x =0 ; x < vol ; ++x)
		fft->rdensity[x] = temp[x] ;
	fft_forward(fft) ;
	
	// Move (q=0) to center of model
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		model[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)] = fft->fdensity[x*size*size + y*size + z] / vol ;
	fprintf(stderr, "Parsed %s: ", fname) ;
	fprintf(stderr, "Center of mass = (%.2f, %.2f, %.2f)\n", cm[0], cm[1], cm[2]) ;
	
	free(temp) ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, bin, num_bins, vox, *numvox ;
	double d_min, r, binsize ;
	float *norm1, *norm2 ;
	float complex *model1, *model2, *fsc ;
	FILE *fp ;
	char fname[999] ;
	float cm1[3] = {0.}, cm2[3] = {0.} ;
	struct fft_data fft ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <model1> <model2> <res_at_edge>\n", argv[0]) ;
		fprintf(stderr, "Optional: <fsc_fname>\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	d_min = atof(argv[3]) ;
	if (argc > 4) {
		strcpy(fname, argv[4]) ;
	}
	else {
		sprintf(fname, "%s", extract_fname(argv[1])) ;
		strtok(fname, "_.") ;
		int num1 = atoi(strtok(NULL, "_.")) ;
		sprintf(fname, "%s", extract_fname(argv[2])) ;
		strtok(fname, "_.") ;
		int num2 = atoi(strtok(NULL, "_.")) ;
		sprintf(fname, "fsc-%d-%d.dat", num1, num2) ;
	}
	
	vol = size*size*size ;
	c = size / 2 ;
	fprintf(stderr, "Model size = %ld\n", size) ;
	num_bins = 150 ;
	binsize = ((double) c) / num_bins ;

	// Allocate memory
	fft_init(&fft, size, omp_get_max_threads()) ;
	fft_create_plans(&fft) ;
	
	model1 = malloc(vol * sizeof(float complex)) ;
	model2 = malloc(vol * sizeof(float complex)) ;
	fsc = calloc(num_bins, sizeof(float complex)) ;
	norm1 = calloc(num_bins, sizeof(float)) ;
	norm2 = calloc(num_bins, sizeof(float)) ;
	numvox = calloc(num_bins, sizeof(long)) ;
	
	process_model(argv[1], model1, &fft, cm1) ;
	process_model(argv[2], model2, &fft, cm2) ;

	// Calculate FSC while correcting for COM difference
	cm2[0] -= cm1[0] ;
	cm2[1] -= cm1[1] ;
	cm2[2] -= cm1[2] ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		r = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		bin = (int) (r / binsize) ;
		if (bin > num_bins)
			continue ;
		
		vox = x*size*size + y*size + z ;
		fsc[bin] += model1[vox] * conjf(model2[vox]) * cexpf(-I*2.*M_PI*((x-c)*cm2[0] + (y-c)*cm2[1] + (z-c)*cm2[2])/size) ;
		norm1[bin] += powf(cabsf(model1[vox]), 2.f) ;
		norm2[bin] += powf(cabsf(model2[vox]), 2.f) ;
		numvox[bin]++ ;
	}
	
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (norm1[bin] * norm2[bin] > 0.)
		fsc[bin] /= sqrtf(norm1[bin] * norm2[bin]) ;

	// Write to file
	fprintf(stderr, "Writing to %s\n", fname) ;
	fp = fopen(fname, "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		fprintf(fp, "%.3f\t%.3f\t%.6f\t%.8ld\n", (bin+1)/d_min/num_bins, num_bins*d_min/(bin+1), cabsf(fsc[bin]), numvox[bin]) ;
	fclose(fp) ;

	// Free memory
	free(model1) ;
	free(model2) ;
	fft_free(&fft) ;
	free(fsc) ;
	free(norm1) ;
	free(norm2) ;
	free(numvox) ;
	
	return 0 ;
}

