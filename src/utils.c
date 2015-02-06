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
	long ch = hsize/2+1, ck = ksize/2+1, cl = lsize/2+1, maxc ;
	float *contrast = calloc(num_bins, sizeof(float)) ;
	long *bin_count = calloc(num_bins, sizeof(long)) ;
	FILE *fp ;
	
	// Continuous part
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = model[x] ;
	
	fftwf_execute(forward_cont) ;
	
	symmetrize_intens(fdensity, exp_intens, 1) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = (x + center1) % size  - center1 ;
		dy = (y + center1) % size  - center1 ;
		dz = (z + center1) % size  - center1 ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / center1 * num_bins + 0.5 ;
		
		if (bin < num_bins && obs_mag[x*size*size + y*size + z] > 0.) {
//			contrast[bin] += cabs(fdensity[x*size*size + y*size + z]) / 
			contrast[bin] += sqrt(exp_intens[x*size*size + y*size + z]) / 
			                 obs_mag[x*size*size + y*size + z] ;
			bin_count[bin]++ ;
		}
	}
	
	// Bragg part
	maxc = ch > ck ? ch : ck ;
	maxc = maxc > cl ? maxc : cl ;
	
	fprintf(stderr, "maxc = %ld\n", maxc) ;
	
	for (x = 0 ; x < hsize ; ++x)
	for (y = 0 ; y < ksize ; ++y)
	for (z = 0 ; z < lsize ; ++z)
		rhkl[x*ksize*lsize + y*lsize + z] 
		 = model[(x+hoffset)*size*size + (y+koffset)*size + (z+loffset)] ;
	
	fftwf_execute(forward_hkl) ;
	
	symmetrize_intens(fhkl, exp_hkl, 0) ;
	
	long mbin = 0 ;
	for (x = 0 ; x < hsize ; ++x)
	for (y = 0 ; y < ksize ; ++y)
	for (z = 0 ; z < lsize ; ++z) {
		dx = (x + ch) % hsize  - ch ;
		dy = (y + ck) % ksize  - ck ;
		dz = (z + cl) % lsize  - cl ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / maxc * num_bins + 0.5 ;
		if (bin > mbin)
			mbin = bin ;
		
		if (bin < num_bins && hkl_mag[x*size*size + y*size + z] > 0.) {
//			contrast[bin] += cabs(fdensity[x*size*size + y*size + z]) / 
			contrast[bin] += sqrt(exp_hkl[x*ksize*lsize + y*lsize + z]) / 
			                 hkl_mag[x*ksize*lsize + y*lsize + z] ;
			bin_count[bin]++ ;
		}
	}
	
	fprintf(stderr, "mbin = %ld\n", mbin) ;
	
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

// Symmetrize intensity according to P(2_1)(2_1)(2_1) space group
// (size,hsize,ksize,lsize)
void symmetrize_intens(fftwf_complex *in, float *out, int flag) {
	long hs, ks, ls, kc, lc ;
	
	if (flag == 1) {
		hs = size ;
		ks = size ;
		ls = size ;
		kc = ks / 2 ;
		lc = ls / 2 ;
	}
	else if (flag == 0) {
		hs = hsize ;
		ks = ksize ;
		ls = lsize ;
		kc = ks / 2 ;
		lc = ls / 2 ;
	}
	
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
		}
	}
}
