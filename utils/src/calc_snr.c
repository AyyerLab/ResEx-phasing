#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <omp.h>
#include "../../src/utils.h"

long size ;

int init_radavg(int *rad, long *rad_count, double binsize) {
	long x, y, z, c = size/2 ;
	int bin, max_bin = 0 ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		bin = floor(sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) / binsize) ;
		rad[x*size*size + y*size + z] = bin ;
		rad_count[bin]++ ;
		if (bin > max_bin)
			max_bin = bin ;
	}
	
	for (bin = 0 ; bin < max_bin ; ++bin)
	if (rad_count[bin] == 0)
		rad_count[bin] = 1 ;
	
	return max_bin ;
}

void calc_radavg(float *model, int *rad, long *count, float *avg) {
	long x, vol = size*size*size ;
	
	for (x = 0 ; x < vol ; ++x)
		avg[rad[x]] += model[x] ;
	
	for (x = 0 ; x < size; ++x)
		avg[x] /= count[x] ;
}

void calc_radial_properties(float *model, int *rad, long *count, float *avg, float *min, float *max) {
	long x, vol = size*size*size ;
	
	for (x = 0 ; x < vol ; ++x) {
		float val = model[x] ;
		int bin = rad[x] ;
		avg[bin] += val ;
		if (val < min[bin])
			min[bin] = val ;
		if (val > max[bin])
			max[bin] = val ;
	}
	
	for (x = 0 ; x < size ; ++x)
		avg[x] /= count[x] ;
}

void calc_histograms(float *model, int *rad, long *count, float *min, float *max, int num_rbins, int num_hbins, double *hist) {
	long x, vol = size*size*size ;
	int rbin, hbin ;
	float val ;
	
	for (x = 0 ; x < vol ; ++x) {
		rbin = rad[x] ;
		if (min[rbin] == max[rbin])
			continue ;
		val = model[x] ;
		hbin = (int) floor((val - min[rbin]) / (max[rbin] - min[rbin]) * num_hbins) ;
		hist[rbin*num_hbins + hbin]++ ;
	}
	
	for (rbin = 0 ; rbin < num_rbins ; ++rbin)
	for (hbin = 0 ; hbin < num_hbins ; ++hbin)
		hist[rbin*num_hbins + hbin] /= (double) count[rbin] ;
}

int main(int argc, char *argv[]) {
	long vol, i, j ;
	float *intens, *sigma, *kl_div ;
	float *radavg, *radsqavg, *radmax, *radmin ;
	int num_bins, *radius, num_histbins = 100 ;
	long *radcount ;
	double binsize, *histograms, p_normal ;
	char fname[999] ;
	FILE *fp ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <intens_fname> <binsize>\n", argv[0]) ;
		fprintf(stderr, "Optional: <output_fname>\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	binsize = atof(argv[2]) ;
	if (argc > 3)
		strcpy(fname, argv[3]) ;
	else
		sprintf(fname, "%s-snr.dat", remove_ext(argv[1])) ;
	vol = size*size*size ;
	
	intens = malloc(vol * sizeof(float)) ;
	radius = malloc(vol * sizeof(int)) ;
	radcount = calloc(size, sizeof(long)) ;
	radavg = calloc(size, sizeof(float)) ;
	radmax = calloc(size, sizeof(float)) ;
	radmin = calloc(size, sizeof(float)) ;
	radsqavg = calloc(size, sizeof(float)) ;
	sigma = calloc(size, sizeof(float)) ;
	kl_div = calloc(size, sizeof(float)) ;
	
	// Read in model
	fp = fopen(argv[1], "rb") ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Parsed intens\n") ;
	
	for (i = 0 ; i < vol ; ++i)
	if (intens[i] < -500.)
		intens[i] = 0. ;
	
	// Calculate radial bins and occupancies
	num_bins = init_radavg(radius, radcount, binsize) ;
	fprintf(stderr, "Calculated radial bins\n") ;
	
	for (i = 0 ; i < num_bins ; ++i) {
		radmax[i] = -FLT_MAX ;
		radmin[i] = FLT_MAX ;
	}
	
	// Calculate radial average, maximum and minimum
	calc_radial_properties(intens, radius, radcount, radavg, radmin, radmax) ;
	
	// Calculate histogram in each radial bin
	histograms = calloc(num_bins * num_histbins, sizeof(double)) ;
	calc_histograms(intens, radius, radcount, radmin, radmax, num_bins, num_histbins, histograms) ;
	fprintf(stderr, "Calculated histograms in each radial bin\n") ;
	
	// Calculate squared radial average and then sigma
	for (i = 0 ; i < vol ; ++i)
		intens[i] = intens[i]*intens[i] ;
	
	calc_radavg(intens, radius, radcount, radsqavg) ;
	
	for (i = 0 ; i < num_bins; ++i)
		sigma[i] = sqrtf(radsqavg[i] - radavg[i]*radavg[i]) ;
	
	// Calculate KL divergence between Gaussian and measured histogram
	for (i = 0 ; i < num_bins ; ++i) {
		double norm_factor = 1. / sqrt(2. * M_PI) / sigma[i] ;
		double hist_binsize = (radmax[i] - radmin[i]) / num_histbins ;
		
		for (j = 0 ; j < num_histbins ; ++j) {
			p_normal = norm_factor * exp(-pow(radmin[i] + j*hist_binsize + 0.5*hist_binsize - radavg[i], 2.) / (2. * sigma[i] * sigma[i])) ;
			if (p_normal > 1.e-6 && histograms[i*num_histbins + j] > 0.) {
				kl_div[i] += p_normal * log(p_normal / histograms[i*num_histbins + j]) ;
				kl_div[i] += histograms[i*num_histbins + j] * log(histograms[i*num_histbins + j] / p_normal) ;
			}
		}
	}
	
/*	int num_sym = 4 ;
	for (i = 0 ; i < num_bins; ++i) {
		sigma[i] = sqrtf(radsqavg[i] - radavg[i]*radavg[i]) ;
		if (isnan(sigma[i]))
			sigma[i] = 0.f ;
		
		kl_div[i] = num_sym * log(num_sym / radavg[i])
		          + log(sqrt(2.*M_PI) * sigma[i] / gsl_sf_fact(num_sym-1))
		          + (num_sym-1)*(gsl_sf_psi_int(num_sym) - log(num_sym/radavg[i]))
				  - num_sym
				  + radavg[i]*radavg[i] / (2. * sigma[i]*sigma[i]) ;
		if (isnan(kl_div[i]))
			kl_div[i] = 0.f ;
	}
*/	
	fprintf(stderr, "Writing output to %s\n", fname) ;
	fp = fopen(fname, "w") ;
	fprintf(fp, "Radius Min    Mean   Max    Std    D_KL\n") ;
	for (i = 0 ; i < num_bins; ++i)
		fprintf(fp, "%6.2f %+.6e %+.6e %.6e %.6e %+.6e\n", binsize*i, radmin[i], radavg[i], radmax[i], sigma[i], kl_div[i]) ;
	fclose(fp) ;
	
	free(intens) ;
	free(radius) ;
	free(radmax) ;
	free(radmin) ;
	free(radcount) ;
	free(radavg) ;
	free(radsqavg) ;
	free(sigma) ;
	free(kl_div) ;
	free(histograms) ;
	
	return 0 ;
}
