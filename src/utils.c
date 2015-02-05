#include "brcont.h"

void init_model(float *model) {
	long i ;
	struct timeval t1 ;
	const gsl_rng_type *T ;
	gsl_rng *r ;
	
	gsl_rng_env_setup() ;
	T = gsl_rng_default ;
	gettimeofday(&t1, NULL) ;
	r = gsl_rng_alloc(T) ;
	gsl_rng_set(r, t1.tv_sec + t1.tv_usec) ;
	
	memset(model, 0, vol*sizeof(float)) ;
	
	for (i = 0 ; i < num_supp ; ++i)
		model[support[i]] = gsl_rng_uniform(r) ;
	
	gsl_rng_free(r) ;
}

void average_model(float *current, float *sum) {
	long i ;
	
	for (i = 0 ; i < vol ; ++i)
		sum[i] += current[i] ;
}

void gen_prtf(float *model) {
	long x, y, z, bin, num_bins = 50 ;
	long dx, dy, dz, center1 = size / 2 + 1 ;
	float scale_factor = sqrt(vol) ;
	float *contrast = calloc(num_bins, sizeof(float)) ;
	long *bin_count = calloc(num_bins, sizeof(long)) ;
	FILE *fp ;
	
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = model[x] ;
	
	fftwf_execute(forward_cont) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = (x + center1) % size  - center1 ;
		dy = (y + center1) % size  - center1 ;
		dz = (z + center1) % size  - center1 ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / center1 * num_bins + 0.5 ;
		
		if (bin < num_bins && obs_mag[x*size*size + y*size + z] > 0.) {
			contrast[bin] += cabs(fdensity[x*size*size + y*size + z]) / 
				(scale_factor * obs_mag[x*size*size + y*size + z]) ;
			bin_count[bin]++ ;
		}
	}
	
	fp = fopen("prtf.dat", "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (bin_count[bin] > 0)
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, contrast[bin]/bin_count[bin]) ;
	else 
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, 0.) ;
	fclose(fp) ;
	
	free(contrast) ;
	free(bin_count) ;
}

