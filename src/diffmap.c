#include "brcont.h"

// Data projection 
// Can set out = in
void proj_data(float *in, float *out) {
	long i ;
	float norm_factor = 1.f / (float) vol ;
	
	// Fourier transform to get structure factors
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = in[i] ;
	
	fftwf_execute(forward) ;
	
	// Replace with known Bragg magnitudes and phases
	for (i = 0 ; i < vol ; ++i)
	if (bragg_calc[i] != FLT_MAX)
		fdensity[i] = bragg_calc[i] ;
//	else                                               // Only when doing Bragg-only reconstruction 
//		fdensity[i] = 0.f ;
	
	// Symmetrize to get intensities to compare
	symmetrize_incoherent(fdensity, exp_mag) ;
	
	// Scale using measured modulus at high resolution
	for (i = 0 ; i < vol ; ++i) {
		if (obs_mag[i] > 0.)
			fdensity[i] *= obs_mag[i] / exp_mag[i] ;
		else if (obs_mag[i] == 0.)
			fdensity[i] = 0. ;
	}
	
	// Inverse Fourier transform
	fftwf_execute(inverse) ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = crealf(rdensity[i]) * norm_factor ;
}

// Support projection
// (vol,num_supp,support)
// Cannot set out = in
void proj_supp(float * restrict in, float * restrict out) {
	long i ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = in[i] * support[i] ;
}

double diffmap(float *x) {
	long i ;
	float diff, change = 0.f ;
//	float alpha = 0.1 ;
//	float beta = 0.7 ;
	
//	proj_data(x, p1) ;                                 // for alpha != 0
	
//	for (i = 0 ; i < vol ; ++i)
//		x[i] = alpha*x[i] + (1.-alpha)*p1[i] ;         // for alpha != 0
	
	proj_supp(x, p1) ;
//	proj_data(x, p2) ;                                 // for beta != 1
	
	for (i = 0 ; i < vol ; ++i) {
		r1[i] = 2. * p1[i] - x[i] ;                    // for beta == 1
//		r1[i] = (1. + 1./beta) * p1[i] - x[i] / beta ; // for beta != 1
//		r2[i] = (1. - 1./beta) * p2[i] + x[i] / beta ; // for beta != 1
	}
	
	proj_data(r1, p2) ;
//	proj_supp(r2, p1) ;                                // for beta != 1
	
	for (i = 0 ; i < vol ; ++i) {
		diff = p2[i] - p1[i] ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

double error_red(float *x) {
	long i ;
	float diff, change = 0.f ;
	
	proj_data(x, p1) ;
	proj_supp(p1, p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = p2[i] - p1[i] ;
		x[i] = p2[i] ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

double modified_hio(float *x) {
	long i ;
	float diff, change = 0.f ;
	float beta = 0.9 ;
	float thresh = 0.1 ;
	
	proj_data(x, p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		r1[i] = (1.f + beta) * p1[i] - x[i] ;
	
	proj_supp(r1, p2) ;
	proj_supp(p1, r1) ;
	
	for (i = 0 ; i < vol ; ++i) {
		if (fabs(p1[i]) > thresh)
			diff = p2[i] - beta*p1[i] ;
		else
			diff = r1[i] - x[i] ;
		
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}
