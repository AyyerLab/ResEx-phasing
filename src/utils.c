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
	long dx, dy, dz, c = size/2, c1 = size/2+1 ;
	float *contrast = calloc(num_bins, sizeof(float)) ;
	long *bin_count = calloc(num_bins, sizeof(long)) ;
	FILE *fp ;
	
	// Continuous part
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = model[x] ;
	
	fftwf_execute(forward) ;
	
	symmetrize_incoherent(fdensity, exp_mag) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		model[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
			= pow(cabsf(fdensity[x*size*size + y*size + z]), 2.) ;
		
		dx = (x + c1) % size - c1 ;
		dy = (y + c1) % size - c1 ;
		dz = (z + c1) % size - c1 ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / c1 * num_bins + 0.5 ;
		
		if (bin < num_bins && obs_mag[x*size*size + y*size + z] > 0.) {
//			contrast[bin] += cabs(fdensity[x*size*size + y*size + z]) / 
			contrast[bin] += exp_mag[x*size*size + y*size + z]
			                 / obs_mag[x*size*size + y*size + z] ;
			bin_count[bin]++ ;
		}
	}
	
	fp = fopen("data/bfrecon.raw", "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	fp = fopen("prtf.dat", "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (bin_count[bin] > 0)
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, contrast[bin]/bin_count[bin]) ;
	else 
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, 1.) ;
	fclose(fp) ;
	
	free(contrast) ;
	free(bin_count) ;
}

// Symmetrize intensity incoherently according to P(2_1)(2_1)(2_1) space group
// (size)
void symmetrize_incoherent(fftwf_complex *in, float *out) {
	long hs, ks, ls, kc, lc ;
	
	hs = size ;
	ks = size ;
	ls = size ;
	kc = ks / 2 ;
	lc = ls / 2 ;
	
	#pragma omp parallel default(shared)
	{
		long x, y, z ;
		float ave_intens ;
		
		#pragma omp for schedule(static,1)
		for (x = 0 ; x < hs ; ++x) {
			for (y = 1 ; y <= kc ; ++y)
			for (z = 1 ; z <= lc ; ++z) {
				ave_intens = 0.25 * (
					powf(cabsf(in[x*ks*ls + y*ls + z]), 2.) +
					powf(cabsf(in[x*ks*ls + y*ls + (ls-z)]), 2.) +
					powf(cabsf(in[x*ks*ls + (ks-y)*ls + z]), 2.) +
					powf(cabsf(in[x*ks*ls + (ks-y)*ls + (ls-z)]), 2.)) ;
				
				out[x*ks*ls + y*ls + z] = ave_intens ;
				out[x*ks*ls + y*ls + (ls-z)] = ave_intens ;
				out[x*ks*ls + (ks-y)*ls + z] = ave_intens ;
				out[x*ks*ls + (ks-y)*ls + (ls-z)] = ave_intens ;
			}
			
			for (z = 1 ; z <= lc ; ++z) {
				ave_intens = 0.5 * (
					powf(cabsf(in[x*ks*ls + z]), 2.) + 
					powf(cabsf(in[x*ks*ls + (ls-z)]), 2.)) ;
				
				out[x*ks*ls + z] = ave_intens ;
				out[x*ks*ls + (ls-z)] = ave_intens ;
			}
			
			for (y = 1 ; y <= kc ; ++y) {
				ave_intens = 0.5 * (
					powf(cabsf(in[x*ks*ls + y*ls]), 2.) + 
					powf(cabsf(in[x*ks*ls + (ks-y)*ls]), 2.)) ;
				
				out[x*ks*ls + y*ls] = ave_intens ;
				out[x*ks*ls + (ks-y)*ls] = ave_intens ;
			}
			
			out[x*ks*ls] = powf(cabsf(in[x*ks*ls]), 2.) ;
			
			for (y = 0 ; y < ks ; ++y)
			for (z = 0 ; z < ls ; ++z)
				out[x*ks*ls + y*ls + z] = sqrtf(out[x*ks*ls + y*ls + z]) ;
		}
	}
}

// Recalculate support
// 'blur' gives width of Gaussian used to convolve with density
// 'threshold' gives cutoff value as a fraction of maximum
void apply_shrinkwrap(float *model, float blur, float threshold) {
	long x, y, z, c = size/2 ;
	float rsq, fblur ;
	uint8_t *supvol = calloc(vol, sizeof(uint8_t)) ;
	
	for (x = 0 ; x < vol ; ++x) {
		rdensity[x] = model[x] ;
	}
	
	// Blur density
	fftwf_execute(forward) ;
	
	fblur = size / (2. * M_PI * blur) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rsq = ((x+c)%size-c)*((x+c)%size-c) + ((y+c)%size-c)*((y+c)%size-c) + ((z+c)%size-c)*((z+c)%size-c) ;
		fdensity[x*size*size + y*size + z] *= expf(-rsq / 2. / fblur / fblur) ;
	}
	
	fftwf_execute(inverse) ;
	
	// Apply threshold
	num_supp = 0 ;
	
	float max = -1.f ;
	for (x = 0 ; x < vol ; ++x)
//	if (crealf(rdensity[x]) > max)
	if (cabsf(rdensity[x]) > max)
//		max = crealf(rdensity[x]) ;
		max = cabsf(rdensity[x]) ; // Absolute support
	
//	threshold *= max ;
	
	for (x = 0 ; x < vol ; ++x)
//	if (crealf(rdensity[x]) > threshold) {
	if (crealf(rdensity[x]) > threshold) { // Absolute support
		support[num_supp++] = x ;
		supvol[x] = 1 ;
	}
	
	char fname[999] ;
	sprintf(fname, "data/shrinkwrap_501_%d.supp", iter) ;
	FILE *fp = fopen(fname, "wb") ;
	fwrite(supvol, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
/*	fp = fopen("data/smoothed.raw", "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
*/	
	free(supvol) ;
}
