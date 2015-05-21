/* Phasing the continuous intensities using Bragg phases and a sequential
	approach. The Bragg projection involves replacing both the Bragg 
	magnitudes and phases for those hkl values where they are known. The 
	continuous projection replaces the magnitudes with the known magnitudes
	and leaves the phases unchanged. The other constraint is the support 
	constraint. The Difference Map update with beta = 1 is used.
*/

#include "brcont.h"

// Bragg projection 
// (hklvol,vol,hsize,ksize,lsize,hoffset,koffset,loffset)
// (rhkl,fhkl,exp_hkl,hkl_mag)
// Can set out = in
void proj_bragg(float *in, float *out) {
	long i ;
	float norm_factor = 1.f / (float) vol ;
	
	// Fourier transform to get structure factors
	for (i = 0 ;i < vol ; ++i)
		rdensity[i] = in[i] ;
	
	fftwf_execute(forward) ;
	
	// Replace with known magnitudes and phases
	for (i = 0 ; i < vol ; ++i)
	if (bragg_calc[i] != FLT_MAX)
		fdensity[i] = bragg_calc[i] ;
//	else /* Only when doing Bragg-only reconstruction */
//		fdensity[i] = 0.f ;
	
	// Inverse Fourier transform
	fftwf_execute(inverse) ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = crealf(rdensity[i]) * norm_factor ;
	
	FILE *fp = fopen("data/bragg_dens.raw", "wb") ;
	fwrite(out, sizeof(float), vol, fp) ;
	fclose(fp) ;
}

// Projection over continuous data 
// (vol,rdensity,fdensity,exp_intens,obs_mag)
// Can set out = in
void proj_cont(float *in, float *out) {
	long i ;
	float norm_factor = 1.f / (float) vol ;
	memset(exp_mag, 0, vol*sizeof(float)) ;
	
	// Fourier transform to get structure factors
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = in[i] ;
	
	fftwf_execute(forward) ;
	
	// Symmetrize to get intensities to compare
	symmetrize_incoherent(fdensity, exp_mag) ;
	
	// Scale using measured modulus at high resolution
	for (i = 0 ; i < vol ; ++i) {
		if (obs_mag[i] > 0.)
			fdensity[i] *= obs_mag[i] / exp_mag[i] ;
		else if (obs_mag[i] == 0.)
			fdensity[i] = 0. ;
/*		else if (obs_mag[i] != -1.f) {
			if (exp_mag[i] > 100.)
				fdensity[i] *= 100. / exp_mag[i] ;
		}
*/	}
	
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
		if (in[pixel] > 0.)
			out[pixel] = in[pixel] ;
	}
}

// Data constraint (Bragg + Continuous)
// Can set out = in
void proj_data(float *in, float *out) {
	long i ;
	
	proj_bragg(in, out) ;
	for (i = 0 ; i < vol ; ++i)
		in[i] = out[i] ;
	proj_cont(in, out) ;
}

double diffmap(float *x) {
	long i ;
	float diff, change = 0. ;
	
	proj_data(x, p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		r1[i] = 2. * p1[i] - x[i] ;
	
	proj_supp(r1, p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = p2[i] - p1[i] ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}
