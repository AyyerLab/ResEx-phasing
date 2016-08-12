#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <omp.h>

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, bin, num_bins, vox, *numvox ;
	double d_min, r, binsize ;
	float *temp, *norm1, *norm2 ;
	float complex *model1, *model2, *fsc ;
	fftwf_complex *rdensity, *fdensity ;
	fftwf_plan forward ;
	FILE *fp ;
	char fname[999] ;
	double cm1[3], cm1_denr ;
	double cm2[3], cm2_denr ;
	
	if (argc < 5) {
		fprintf(stderr, "Format: %s <model1> <model2> <size> <res_at_edge>\n", argv[0]) ;
		fprintf(stderr, "Optional: <fsc_fname>\n") ;
		return 1 ;
	}
	size = atoi(argv[3]) ;
	d_min = atof(argv[4]) ;
	vol = size*size*size ;
	c = size / 2 ;
	fprintf(stderr, "Model size = %ld\n", size) ;
	
	// Allocate memory
	temp = malloc(vol * sizeof(float)) ;
	rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	model1 = malloc(vol * sizeof(float complex)) ;
	model2 = malloc(vol * sizeof(float complex)) ;
	
	num_bins = 150 ;
	binsize = ((double) c) / num_bins ;
	
	fsc = calloc(num_bins, sizeof(float complex)) ;
	norm1 = calloc(num_bins, sizeof(float)) ;
	norm2 = calloc(num_bins, sizeof(float)) ;
	numvox = calloc(num_bins, sizeof(long)) ;
	
	// Parse fftwf plan
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(omp_get_max_threads()) ;
	
	forward = fftwf_plan_dft_3d(size, size, size, rdensity, fdensity, FFTW_FORWARD, FFTW_ESTIMATE) ;

	// Parse first model
	fp = fopen(argv[1], "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	// Calculate center of mass
	cm1[0] = 0 ; cm1[1] = 0 ; cm1[2] = 0 ; cm1_denr = 0 ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		cm1[0] += x*temp[x*size*size + y*size + z] ;
		cm1[1] += y*temp[x*size*size + y*size + z] ;
		cm1[2] += z*temp[x*size*size + y*size + z] ;
		cm1_denr += temp[x*size*size + y*size + z] ;
	}
	cm1[0] /= cm1_denr ; cm1[1] /= cm1_denr ; cm1[2] /= cm1_denr ;
	// Calculate Fourier transform
	for (x =0 ; x < vol ; ++x)
		rdensity[x] = temp[x] ;
	fftwf_execute(forward) ;
	// Move (q=0) to center of model
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		model1[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)] = fdensity[x*size*size + y*size + z] / vol ;
	fprintf(stderr, "Parsed %s: ", argv[1]) ;
	fprintf(stderr, "center of mass = (%.2f, %.2f, %.2f)\n", cm1[0], cm1[1], cm1[2]) ;

	// Parse second model
	fp = fopen(argv[2], "rb") ;
	fread(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	// Calculate center of mass
	cm2[0] = 0 ; cm2[1] = 0 ; cm2[2] = 0 ; cm2_denr = 0 ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		cm2[0] += x*temp[x*size*size + y*size + z] ;
		cm2[1] += y*temp[x*size*size + y*size + z] ;
		cm2[2] += z*temp[x*size*size + y*size + z] ;
		cm2_denr += temp[x*size*size + y*size + z] ;
	}
	cm2[0] /= cm2_denr ; cm2[1] /= cm2_denr ; cm2[2] /= cm2_denr ;
	// Calculate Fourier transform
	for (x =0 ; x < vol ; ++x)
		rdensity[x] = temp[x] ;
	fftwf_execute(forward) ;
	// Move (q=0) to center of model
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		model2[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)] = fdensity[x*size*size + y*size + z] ;
	fprintf(stderr, "Parsed %s: ", argv[2]) ;
	fprintf(stderr, "center of mass = (%.2f, %.2f, %.2f)\n", cm2[0], cm2[1], cm2[2]) ;
	cm2[0] -= cm1[0] ;
	cm2[1] -= cm1[1] ;
	cm2[2] -= cm1[2] ;

	// Calculate and save phase difference while correcting for COM difference
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		temp[x*size*size + y*size + z] = carg(model1[x*size*size + y*size + z]) - carg(model2[x*size*size + y*size + z]) - 2.*M_PI*((x-c)*cm2[0] + (y-c)*cm2[1] + (z-c)*cm2[2])/size ;
	fp = fopen("data/delta_phase.raw", "wb") ;
	fwrite(temp, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate FSC while correcting for COM difference
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
	if (argc > 5) {
		strcpy(fname, argv[5]) ;
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
	
	fprintf(stderr, "Writing to %s\n", fname) ;
	fp = fopen(fname, "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		fprintf(fp, "%.3f\t%.3f\t%.6f\t%.8ld\n", (bin+1)/d_min/num_bins, num_bins*d_min/(bin+1), cabsf(fsc[bin]), numvox[bin]) ;
	fclose(fp) ;
	
	// Free memory
	free(temp) ;
	free(model1) ;
	free(model2) ;
	fftwf_free(rdensity) ;
	fftwf_free(fdensity) ;
	
	return 0 ;
}
