#include "brcont.h"

// Data projection 
// Can set out = in
void proj_data(float *in, float *out) {
	long i, x, y, z, s = size, c = s/2 ;
	float norm_factor = 1.f / (float) vol, n_cells ;
	float mean_incoh = 0, mean_coh = 0, mean_incoh_b = 0 ;
	float mean_obs = 0, mean_obs_b = 0 ;
	long count_nb = 0, count_b = 0 ;
	FILE *fp ;
	
	// Fourier transform to get structure factors
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = in[i] ;
	
	fftwf_execute(forward) ;
	
	// Symmetrize to get intensities to compare
//	symmetrize_incoherent(fdensity, incoh_mag) ;
//	symmetrize_coherent(fdensity, coh_mag) ;
	
/*	FILE *fp = fopen("data/incoh_mag.raw", "wb") ;
	fwrite(incoh_mag, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	fp = fopen("data/coh_mag.raw", "wb") ;
	fwrite(coh_mag, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	fp = fopen("data/obs_mag.raw", "wb") ;
	fwrite(obs_mag, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	fp = fopen("data/fdensity.cpx", "wb") ;
	fwrite(fdensity, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
*/	
	// Scale non-Bragg voxels and calculate n_cells
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z) {
		i = x*s*s + y*s + z ;
		
		if (mask[i] == NON_BRAGG) {
			if (obs_mag[i] > 0.) {
//				mean_incoh += incoh_mag[i] * incoh_mag[i] ;
//				mean_obs += obs_mag[i] * obs_mag[i] ;
//				count_nb++ ;
				
//				fdensity[((x+c)%s)*s*s + ((y+c)%s)*s + ((z+c)%s)] *= obs_mag[i] / incoh_mag[i] ;
				fdensity[((x+c)%s)*s*s + ((y+c)%s)*s + ((z+c)%s)] *= obs_mag[i] / cabsf(fdensity[((x+c)%s)*s*s + ((y+c)%s)*s + ((z+c)%s)]) ;
			}
			else if (obs_mag[i] == 0.) {
				//count_nb++ ;
				fdensity[((x+c)%s)*s*s + ((y+c)%s)*s + ((z+c)%s)] = 0. ;
			}
		}
/*		else if (mask[i] == BRAGG) {
			if (obs_mag[i] > 0. && coh_mag[i] > 0.) {
				mean_incoh_b += incoh_mag[i] * incoh_mag[i] ;
				mean_coh += coh_mag[i] * coh_mag[i] ;
				mean_obs_b += obs_mag[i] * obs_mag[i] ;
				count_b++ ;
			}
		}
*/	}
	
/*	mean_incoh /= count_nb ;
	mean_obs /= count_nb ;
	mean_coh /= count_b ;
	mean_incoh_b /= count_b ;
	mean_obs_b /= count_b ;
	n_cells = (mean_incoh*mean_obs_b/mean_obs - mean_incoh_b) / mean_coh ;
	
	//fprintf(stderr, "means: %.6e, %.6e, %.6e\n", mean_incoh, mean_incoh_b, mean_coh) ;
	//fprintf(stderr, "mean_obs: %.6e, %.6e\n", mean_obs, mean_obs_b) ;
	fprintf(stderr, "n_cells = %f\n", n_cells) ;
*/	n_cells = 1000. ;
	
/*	// Scale Bragg voxels
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z) {
		i = x*s*s + y*s + z ;
		
		if (mask[i] == BRAGG && obs_mag[i] > 0.)
			fdensity[((x+c)%s)*s*s + ((y+c)%s)*s + ((z+c)%s)] 
			  *= obs_mag[i] / sqrtf(n_cells*coh_mag[i]*coh_mag[i] + incoh_mag[i]*incoh_mag[i]) ;
	}
*/	
	// Inverse Fourier transform
	fftwf_execute(inverse) ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = crealf(rdensity[i]) * norm_factor ;
	
	fp = fopen("data/data_proj.raw", "wb") ;
	fwrite(out, sizeof(float), vol, fp) ;
	fclose(fp) ;

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
	float diff, change = 0.f ;
/*	float alpha = 0.95 ;
	
	proj_data(x, p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		x[i] = alpha*x[i] + (1.-alpha)*p1[i] ;
*/	
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
