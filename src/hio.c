#include <brcont.h>

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

double hio(float *x) {
	long i ;
	double change = 0., beta = 1. ;

	proj_data(x, p1) ;
	
	for (i = 0 ; i < vol ; ++i)
	if (supvol[i]) {
		diff = p1[i] - x[i] ;
		x[i] += diff ;
		change += diff*diff ;
	}
	else {
		diff = - beta * p1[x] ;
		x[i] += diff ;
		change += diff*diff ;
	}

	return sqrt(change/vol) ;
}
