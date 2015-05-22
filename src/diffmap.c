#include "brcont.h"

// Data projection 
// Can set out = in
void proj_data(float *in, float *out) {
	long i ;
	float norm_factor = 1.f / (float) vol ;
	
	// Fourier transform to get structure factors
	for (i = 0 ;i < vol ; ++i)
		rdensity[i] = in[i] ;
	
	fftwf_execute(forward) ;
	
	// Replace with known Bragg magnitudes and phases
	for (i = 0 ; i < vol ; ++i)
	if (bragg_calc[i] != FLT_MAX)
		fdensity[i] = bragg_calc[i] ;
//	else // Only when doing Bragg-only reconstruction 
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
void proj_supp(float *in, float *out) {
	long i, pixel ;
	
	memset(out, 0, vol*sizeof(float)) ;
	
	for (i = 0 ; i < num_supp ; ++i) {
		pixel = support[i] ;
//		if (in[pixel] > 0.) // Positivity
			out[pixel] = in[pixel] ;
	}
}

double diffmap(float *x) {
	long i ;
	float diff, change = 0. ;
	
	proj_supp(x, p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		r1[i] = 2. * p1[i] - x[i] ;
	
	proj_data(r1, p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = p2[i] - p1[i] ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}
