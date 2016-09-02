#include "brcont.h"

/* Fourier-space projection 
 * Makes model consistent with data. Modifies Fourier magnitudes and keeps
 * phases unchanged. Applies symmetry depending on the specified point group.
 * Can be applied in-place (out = in)
 */
void proj_fourier(float *in, float *out) {
	long i ;
	float norm_factor = 1.f / (float) vol ;
	
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = in[i] ;
	
	fftwf_execute(forward_plan) ;
	
	symmetrize_incoherent(fdensity, exp_mag) ;
	
	for (i = 0 ; i < vol ; ++i)
	if (bragg_calc[i] != FLT_MAX)
		fdensity[i] = bragg_calc[i] ;
	
	for (i = 0 ; i < vol ; ++i) {
		if (obs_mag[i] > 0.)
			fdensity[i] *= obs_mag[i] / exp_mag[i] ;
		else if (obs_mag[i] == 0.)
			fdensity[i] = 0. ;
		else if (exp_mag[i] > mag_thresh)
			fdensity[i] *= mag_thresh / cabsf(fdensity[i]) ;
	}
	
	fftwf_execute(inverse_plan) ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = crealf(rdensity[i]) * norm_factor ;
}

/* Direct-space projection
 * Applies the finite support constraint in direct/real-space
 * In addition one can apply the following additional constraints by
 * setting the following global variables.
 * 	do_positivity - Set negative values inside support to zero
 * 	do_histogram - Project values inside support such that they match a 
 * 	               target histogram.
 * Can be applied in-place (out = in)
 */
void proj_direct(float * restrict in, float * restrict out) {
	if (do_histogram) {
		match_histogram(in, out) ;
	}
	else {
		long i ;
		for (i = 0 ; i < vol ; ++i)
			out[i] = in[i] * support[i] ;
		
		if (do_positivity) {
			for (i = 0 ; i < vol ; ++i)
			if (out[i] < 0.)
				out[i] = 0. ;
		}
	}
}

/* Difference Map algorithm
 * Update rule (LaTeX syntax)
 * x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right)P_F(x_n) - \frac{x_n}{\beta}\right] - P_F\left[\left(1-\frac{1}{\beta}\right)P_D(x_n) + \frac{x_n}{\beta}\right]\right\}
 * Same as all other algorithms for beta = 1
 */
double DM_algorithm(float *x) {
	long i ;
	float diff, change = 0.f ;
	
	proj_fourier(x, algorithm_p1) ;
	if (algorithm_beta != 1.)
		proj_direct(x, algorithm_p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		algorithm_r1[i] = (1. + 1./algorithm_beta) * algorithm_p1[i] - x[i] / algorithm_beta ;
		if (algorithm_beta != 1.)
			algorithm_r2[i] = (1. - 1./algorithm_beta) * algorithm_p2[i] + x[i] / algorithm_beta ;
	}
	
	proj_direct(algorithm_r1, algorithm_p2) ;
	if (algorithm_beta != 1.)
		proj_fourier(algorithm_r2, algorithm_p1) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = algorithm_beta * (algorithm_p2[i] - algorithm_p1[i]) ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

/* Modified Difference Map algorithm
 * Update rule (LaTeX syntax)
 * x_n' = \beta x_n + (1-\beta) P_F(x_n)
 * x_{n+1} = x_n' + P_F\left[2 P_D(x_n') - x_n'\right] - P_D(x_n')
 * Same as all other algorithms for beta = 1
 */
double mod_DM_algorithm(float *x) {
	long i ;
	float diff, change = 0.f ;
	
	if (algorithm_beta != 1.) {
		proj_fourier(x, algorithm_p1) ;
	
		for (i = 0 ; i < vol ; ++i)
			x[i] = algorithm_beta*x[i] + (1. - algorithm_beta) * algorithm_p1[i] ;
	}
	
	proj_direct(x, algorithm_p2) ;
	
	for (i = 0 ; i < vol ; ++i)
		algorithm_r1[i] = 2. * algorithm_p2[i] - x[i] ;
	
	proj_fourier(algorithm_r1, algorithm_p1) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = algorithm_p1[i] - algorithm_p2[i] ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

/* RAAR algorithm
 * Update rule (LaTeX syntax)
 * x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n) - x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)
 *
 * If one does not assume P_D is linear,
 * x_{n+1} = \beta \left\{x_n + P_D\left[2 P_F(x_n)\right] + P_D\left[-x_n\right] - P_F(x_n)\right\} + (1-\beta) P_F(x_n)
 * 
 * Same as all other algorithms for beta = 1
 */
double RAAR_algorithm(float *x) {
	long i ;
	float diff, change = 0.f ;
	
	proj_fourier(x, algorithm_p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		algorithm_r1[i] = 2. * algorithm_p1[i] ;
	proj_direct(algorithm_r1, algorithm_r2) ;
	
	for (i = 0 ; i < vol ; ++i)
		algorithm_r1[i] = - algorithm_p1[i] ;
	proj_direct(algorithm_r1, algorithm_p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = (algorithm_beta - 1.) * x[i] + algorithm_beta * (algorithm_r2[i] + algorithm_p2[i]) + (1. - 2. * algorithm_beta) * algorithm_p1[i] ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

/* HIO algorithm
 * Update rule (LaTeX syntax)
 * x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right) P_F(x_n) - \frac{x_n}{\beta}\right] - P_F(x_n)\right\}
 * Same as all other algorithms for beta = 1
 */
double HIO_algorithm(float *x) {
	long i ;
	float diff, change = 0.f ;
	
	proj_fourier(x, algorithm_p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		algorithm_r1[i] = (1. + 1./algorithm_beta) * algorithm_p1[i] - x[i] / algorithm_beta ;
	
	proj_direct(algorithm_r1, algorithm_p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = algorithm_beta * (algorithm_p2[i] - algorithm_p1[i]) ;
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

/* Error Reduction algorithm
 * Update rule (LaTeX style)
 * x_{n+1} = P_D[P_F(x_n)]
 * Obviously different from others. Use only in averaging phase.
 */
double ER_algorithm(float *x) {
	long i ;
	float diff, change = 0.f ;
	
	proj_fourier(x, algorithm_p1) ;
	proj_direct(algorithm_p1, algorithm_p2) ;
	
	for (i = 0 ; i < vol ; ++i) {
		diff = algorithm_p2[i] - algorithm_p1[i] ;
		x[i] = algorithm_p2[i] ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}

double modified_hio(float *x) {
	long i ;
	float diff, change = 0.f ;
	float thresh = 0.1 ;
	
	proj_fourier(x, algorithm_p1) ;
	
	for (i = 0 ; i < vol ; ++i)
		algorithm_r1[i] = (1.f + algorithm_beta) * algorithm_p1[i] - x[i] ;
	
	proj_direct(algorithm_r1, algorithm_p2) ;
	proj_direct(algorithm_p1, algorithm_r1) ;
	
	for (i = 0 ; i < vol ; ++i) {
		if (fabs(algorithm_p1[i]) > thresh)
			diff = algorithm_p2[i] - algorithm_beta*algorithm_p1[i] ;
		else
			diff = algorithm_r1[i] - x[i] ;
		
		x[i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / vol) ;
}
