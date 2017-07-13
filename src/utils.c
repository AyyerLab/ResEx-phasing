#include "brcont.h"

/* Initial guess for model and background
 * Model is white noise inside support volume
 * Background is white noise
*/
void init_model(float *model, int random_model, int init_bg) {
	long i ;
	float val = sqrtf(vol) ;
	struct timeval t1 ;
	const gsl_rng_type *T ;
	gsl_rng *r ;
	
	gsl_rng_env_setup() ;
	T = gsl_rng_default ;
	gettimeofday(&t1, NULL) ;
	r = gsl_rng_alloc(T) ;
	gsl_rng_set(r, t1.tv_sec + t1.tv_usec) ;
	
	if (random_model)
		memset(model, 0, vol*sizeof(float)) ;
	if (init_bg)
		memset(&(model[vol]), 0, vol*sizeof(float)) ;
	
	if (random_model || init_bg) {
		for (i = 0 ; i < vol ; ++i)
		if (random_model && support[i])
			model[i] = gsl_rng_uniform(r) ;
		
		if (do_bg_fitting) {
			for (i = 0 ; i < vol ; ++i)
			if (init_bg && obs_mag[i] > 0.)
				model[vol+i] = val ;
		}
	}
	
	gsl_rng_free(r) ;
}

void average_model(float *current, float *sum) {
	long i, num_vox = vol ;
	
	if (do_bg_fitting)
		num_vox *= 2 ;
	
	for (i = 0 ; i < num_vox ; ++i)
		sum[i] += current[i] ;
}

/* Calculate PRTF and save estimated intensities due to
*/
void gen_prtf(float *model) {
	long x, y, z, bin, num_bins = 50 ;
	long dx, dy, dz, c = size/2, c1 = size/2+1 ;
	float obs_val, *contrast = calloc(num_bins, sizeof(float)) ;
	long *bin_count = calloc(num_bins, sizeof(long)) ;
	FILE *fp ;
	char fname[999] ;
	
	// Continuous part
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = model[x] ;
	
	fftwf_execute(forward_plan) ;
	
	if (do_bg_fitting)
		symmetrize_incoherent(fdensity, exp_mag, &(model[vol])) ;
	else
		symmetrize_incoherent(fdensity, exp_mag, NULL) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		model[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
			= pow(cabsf(fdensity[x*size*size + y*size + z]), 2.) ;
		
		dx = (x + c1) % size - c1 ;
		dy = (y + c1) % size - c1 ;
		dz = (z + c1) % size - c1 ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / c1 * num_bins + 0.5 ;
		obs_val = obs_mag[x*size*size + y*size + z] ;
		
		if (bin < num_bins && obs_val > 0.) {
			contrast[bin] += exp_mag[x*size*size + y*size + z]
			                 / obs_mag[x*size*size + y*size + z] ;
			bin_count[bin]++ ;
		}
		else if (bin < num_bins && obs_val == 0.)
			bin_count[bin]++ ;
	}
	
	sprintf(fname, "%s-expmag.raw", output_prefix) ;
	dump_slices(exp_mag, fname, 1) ;
	
	sprintf(fname, "%s-frecon.raw", output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	sprintf(fname, "%s-prtf.dat", output_prefix) ;
	fp = fopen(fname, "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (bin_count[bin] > 0)
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, contrast[bin]/bin_count[bin]) ;
	else 
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, 1.) ;
	fclose(fp) ;
	
	free(contrast) ;
	free(bin_count) ;
}

/* Symmetrize intensity incoherently according to 222 point group
 * The array is assumed to have q=0 at (0,0,0) instead of in the center of the array
 * Global variables: (size, point_group)
*/
void symmetrize_incoherent(fftwf_complex *in, float *out, float *bg) {
	long hs, ks, ls, kc, lc ;
	
	hs = size ;
	ks = size ;
	ls = size ;
	kc = ks / 2 ;
	lc = ls / 2 ;
	
	if (strcmp(point_group, "222") == 0) {
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
					
					if (bg != NULL) {
						ave_intens = 0.25 * (
							powf(bg[x*ks*ls + y*ls + z], 2.) +
							powf(bg[x*ks*ls + y*ls + (ls-z)], 2.) +
							powf(bg[x*ks*ls + (ks-y)*ls + z], 2.) +
							powf(bg[x*ks*ls + (ks-y)*ls + (ls-z)], 2.)) ;
						
						bg[x*ks*ls + y*ls + z] = ave_intens ;
						bg[x*ks*ls + y*ls + (ls-z)] = ave_intens ;
						bg[x*ks*ls + (ks-y)*ls + z] = ave_intens ;
						bg[x*ks*ls + (ks-y)*ls + (ls-z)] = ave_intens ;
					}
				}
				
				for (z = 1 ; z <= lc ; ++z) {
					ave_intens = 0.5 * (
						powf(cabsf(in[x*ks*ls + z]), 2.) + 
						powf(cabsf(in[x*ks*ls + (ls-z)]), 2.)) ;
					
					out[x*ks*ls + z] = ave_intens ;
					out[x*ks*ls + (ls-z)] = ave_intens ;
					
					if (bg != NULL) {
						ave_intens = 0.5 * (
							powf(bg[x*ks*ls + z], 2.) + 
							powf(bg[x*ks*ls + (ls-z)], 2.)) ;
						
						bg[x*ks*ls + z] = ave_intens ;
						bg[x*ks*ls + (ls-z)] = ave_intens ;
					}
				}
				
				for (y = 1 ; y <= kc ; ++y) {
					ave_intens = 0.5 * (
						powf(cabsf(in[x*ks*ls + y*ls]), 2.) + 
						powf(cabsf(in[x*ks*ls + (ks-y)*ls]), 2.)) ;
					
					out[x*ks*ls + y*ls] = ave_intens ;
					out[x*ks*ls + (ks-y)*ls] = ave_intens ;
					
					if (bg != NULL) {
						ave_intens = 0.5 * (
							powf(bg[x*ks*ls + y*ls], 2.) + 
							powf(bg[x*ks*ls + (ks-y)*ls], 2.)) ;
						
						bg[x*ks*ls + y*ls] = ave_intens ;
						bg[x*ks*ls + (ks-y)*ls] = ave_intens ;
					}
				}
				
				out[x*ks*ls] = powf(cabsf(in[x*ks*ls]), 2.) ;
				if (bg != NULL)
					bg[x*ks*ls] = powf(bg[x*ks*ls], 2.) ;
				
				for (y = 0 ; y < ks ; ++y)
				for (z = 0 ; z < ls ; ++z) {
					if (bg != NULL)
						out[x*ks*ls + y*ls + z] = sqrtf(out[x*ks*ls + y*ls + z] + bg[x*ks*ls + y*ls + z]) ;
					else
						out[x*ks*ls + y*ls + z] = sqrtf(out[x*ks*ls + y*ls + z]) ;
				}
			}
		}
	}
	else if (strcmp(point_group, "4") == 0) {
		#pragma omp parallel default(shared)
		{
			long x, y, z ;
			float ave_intens ;
			
			#pragma omp for schedule(static)
			for (x = 0 ; x < size ; ++x) {
				for (y = 1 ; y <= kc ; ++y)
				for (z = 1 ; z <= lc ; ++z) {
					ave_intens = 0.25 * (
						powf(cabsf(in[x*ks*ls + y*ls + z]), 2.f) +
						powf(cabsf(in[x*ks*ls + (ks-z)*ls + y]), 2.f) +
						powf(cabsf(in[x*ks*ls + z*ls + (ls-y)]), 2.f) + 
						powf(cabsf(in[x*ks*ls + (ks-y)*ls + (ls-z)]), 2.f)) ;
					
					out[x*ks*ls + y*ls + z] = ave_intens ;
					out[x*ks*ls + (ks-z)*ls + y] = ave_intens ;
					out[x*ks*ls + z*ls + (ls-y)] = ave_intens ;
					out[x*ks*ls + (ks-y)*ls + (ls-z)] = ave_intens ;
				}
				
				for (y = 1 ; y <= kc ; ++y) {
					ave_intens = 0.25 * (
						powf(cabsf(in[x*ks*ls + y*ls]), 2.) +
						powf(cabsf(in[x*ks*ls + (ks-y)*ls]), 2.) +
						powf(cabsf(in[x*ks*ls + (ks-y)]), 2.) +
						powf(cabsf(in[x*ks*ls + y]), 2.)) ;
					
					out[x*ks*ls + y*ls] = ave_intens ;
					out[x*ks*ls + (ks-y)*ls] = ave_intens ;
					out[x*ks*ls + (ks-y)] = ave_intens ;
					out[x*ks*ls + y] = ave_intens ;
				}
				
				out[x*ks*ls] = powf(cabsf(in[x*ks*ls]), 2.f) ;
				
				for (y = 0 ; y < ks ; ++y)
				for (z = 0 ; z < ls ; ++z)
				if (bg != NULL)
					out[x*ks*ls + y*ls + z] = sqrtf(out[x*ks*ls + y*ls + z] + powf(bg[x*ks*ls + y*ls + z], 2.f)) ;
				else
					out[x*ks*ls + y*ls + z] = sqrtf(out[x*ks*ls + y*ls + z]) ;
			}
		}
	}
	else if (strcmp(point_group, "1") == 0) {
		long x ;
		for (x = 0 ; x < vol ; ++x)
		if (bg != NULL)
			out[x] = sqrtf(powf(cabsf(in[x]), 2.f) + powf(bg[x], 2.f)) ;
		else
			out[x] = cabsf(in[x]) ;
	}
	else {
		fprintf(stderr, "Unrecognized point group: %s\n", point_group) ;
	}
}

/* Recalculate support
 * 'blur' gives width of Gaussian used to convolve with density
 * 'threshold' gives cutoff value as a fraction of maximum
*/
void apply_shrinkwrap(float *model, float blur, float threshold) {
	long x, y, z, c = size/2 ;
	float rsq, fblur ;
	
	for (x = 0 ; x < vol ; ++x) {
		rdensity[x] = model[x] ;
	}
	
	// Blur density
	fftwf_execute(forward_plan) ;
	
	fblur = size / (2. * M_PI * blur) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rsq = ((x+c)%size-c)*((x+c)%size-c) + ((y+c)%size-c)*((y+c)%size-c) + ((z+c)%size-c)*((z+c)%size-c) ;
		fdensity[x*size*size + y*size + z] *= expf(-rsq / 2. / fblur / fblur) ;
	}
	
	fftwf_execute(inverse_plan) ;
	
	// Apply threshold
	float max = -1.f ;
	for (x = 0 ; x < vol ; ++x)
//	if (crealf(rdensity[x]) > max)
	if (cabsf(rdensity[x]) > max)
//		max = crealf(rdensity[x]) ;
		max = cabsf(rdensity[x]) ; // Absolute support
	
//	threshold *= max ;
	
	for (x = 0 ; x < vol ; ++x)
//	if (crealf(rdensity[x]) > threshold)
	if (crealf(rdensity[x]) > threshold) // Absolute support
		support[x] = 1 ;
	
	char fname[999] ;
	sprintf(fname, "%s-shrink", output_prefix) ;
	mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	sprintf(fname, "%s-shrink/%.4d.supp", output_prefix, iter) ;
	FILE *fp = fopen(fname, "wb") ;
	fwrite(support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
/*	fp = fopen("data/smoothed.raw", "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
*/	
}

/* Save orthogonal central slices to file
*/
void dump_slices(float *vol, char *fname, int is_fourier) {
	long x, y, c = size/2 ;
	FILE *fp ;
	float *slices = malloc(3*size*size*sizeof(float)) ;
	
	if (is_fourier) {
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y) {
			slices[x*size + y] = vol[((x+c)%size)*size + ((y+c)%size)] ;
			slices[size*size + x*size + y] = vol[((x+c)%size)*size*size + ((y+c)%size)] ;
			slices[2*size*size + x*size + y] = vol[((x+c)%size)*size*size + ((y+c)%size)*size] ;
		}
	}
	else {
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y) {
			slices[x*size + y] = vol[c*size*size + x*size + y] ;
			slices[size*size + x*size + y] = vol[x*size*size + c*size + y] ;
			slices[2*size*size + x*size + y] = vol[x*size*size + y*size + c] ;
		}
	}
	
	fp = fopen(fname, "wb") ;
	fwrite(slices, sizeof(float), 3*size*size, fp) ;
	fclose(fp) ;
	
	free(slices) ;
}

void dump_support_slices(uint8_t *vol, char *fname) {
	long x, y, c = size/2 ;
	FILE *fp ;
	uint8_t *slices = malloc(3*size*size*sizeof(float)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y) {
		slices[x*size + y] = vol[c*size*size + x*size + y] ;
		slices[size*size + x*size + y] = vol[x*size*size + c*size + y] ;
		slices[2*size*size + x*size + y] = vol[x*size*size + y*size + c] ;
	}
	
	fp = fopen(fname, "wb") ;
	fwrite(slices, sizeof(uint8_t), 3*size*size, fp) ;
	fclose(fp) ;
	
	free(slices) ;
}

/* Radial average initialization
 * Calculate bin for each voxel and bin occupancy
 * Note that q=0 is at (0,0,0) and not in the center of the array
*/
void init_radavg() {
	long x, y, z, c = size/2 ;
	long dx, dy, dz ;
	double dist ;
	
	intrad = malloc(vol * sizeof(int)) ;
	radavg = calloc(size, sizeof(double)) ;
	obs_radavg = calloc(size, sizeof(double)) ;
	radcount = calloc(size, sizeof(double)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = x <= c ? x : size-x ; 
		dy = y <= c ? y : size-y ; 
		dz = z <= c ? z : size-z ; 
		dist = sqrt(dx*dx + dy*dy + dz*dz) ;
		intrad[x*size*size + y*size + z] = (int) dist ;
//		if (obs_mag[x*size*size + y*size + z] > 0.)
			radcount[intrad[x*size*size + y*size + z]] += 1. ;
	}
	
	for (x = 0 ; x < size ; ++x)
	if (radcount[x] == 0.)
		radcount[x] = 1. ;
}

/* Radial average calculation
 * Using previously calculated bins and bin occupancies
 * Note that q=0 is at (0,0,0) and not in the center of the array
*/
void radial_average(float *in, float *out) {
	long i ;
	
	memset(radavg, 0, size*sizeof(float)) ;
	for (i = 0 ; i < vol ; ++i)
		radavg[intrad[i]] += in[i] ;
	
	for (i = 0 ; i < size ; ++i) {
		radavg[i] /= radcount[i] ;
		if (radavg[i] < 0.)
			radavg[i] = 0. ;
	}
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = radavg[intrad[i]] ;
}

int compare_indices_histogram(const void *a, const void *b) {
	int ia = *(int*) a ;
	int ib = *(int*) b ;
	return supp_val[ia] < supp_val[ib] ? -1 : supp_val[ia] > supp_val[ib] ;
}

/* Histogram matching
 * Matches histogram within support volume using inverse_cdf array.
 * Can be done in-place
 */
void match_histogram(float *in, float *out) {
	long i ; 
	
	for (i = 0 ; i < num_supp ; ++i) {
		supp_index[i] = i ;
		supp_val[i] = in[supp_loc[i]] ;
	}
	
	qsort(supp_index, num_supp, sizeof(long), compare_indices_histogram) ;
	
	memset(out, 0, vol*sizeof(float)) ;
	for (i = 0 ; i < num_supp ; ++i)
		out[supp_loc[supp_index[i]]] = inverse_cdf[i] ;
}

/*Get mode of positive values with support volume
 */
float positive_mode(float *model) {
	long i, valbin, maxhist = 0, hist[99] ;
	float bin[99], maxval = 0. ;
	
	for (i = 0 ; i < vol ; ++i)
	if (model[i] > maxval)
		maxval = model[i] ;
	
	for (i = 0 ; i < 99 ; ++i) {
		hist[i] = 0. ;
		bin[i] = maxval * i / 99. ;
	}
	
	for (i = 0 ; i < vol ; ++i)
	if (model[i] > 0.) {
		valbin = model[i] * 99. / maxval ;
		hist[valbin]++ ;
	}
	
	valbin = 0 ;
	for (i = 0 ; i < 99 ; ++i)
	if (hist[i] > maxhist) {
		valbin = i ;
		maxhist = hist[i] ;
	}
	
	fprintf(stderr, "Mode of positive values in volume = %.3e +- %.3e\n", bin[valbin], maxval / 2. / 99.) ;
	
	return 0.1 * bin[valbin] ;
}

int compare_indices_variation(const void *a, const void *b) {
	int ia = *(int*) a ;
	int ib = *(int*) b ;
	return local_variation[ia] > local_variation[ib] ? -1 : local_variation[ia] < local_variation[ib] ;
}

/* Local variation support update
 * Similar to solvent flattening, keeping number of support voxels the same
 */
void variation_support(float *model, uint8_t *supp, long box_rad) {
	long i, num_vox = 0 ;
	
	memset(local_variation, 0, vol*sizeof(float)) ;
	#pragma omp parallel default(shared)
	{
		int np = omp_get_num_threads(), rank = omp_get_thread_num() ;
		long x, y, z, i, j, k, num_p = 0 ;
		float val, avg, rms, weight = 1. / powf(2 * box_rad + 1, 3.f) ;
		
		#pragma omp for schedule(static,1)
		for (x = support_bounds[0]-3*box_rad ; x <= support_bounds[1]+3*box_rad ; ++x)
		for (y = support_bounds[2]-3*box_rad ; y <= support_bounds[3]+3*box_rad ; ++y)
		for (z = support_bounds[4]-3*box_rad ; z <= support_bounds[5]+3*box_rad ; ++z) {
			avg = 0. ;
			rms = 0. ;
			for (i = -box_rad ; i < box_rad ; ++i)
			for (j = -box_rad ; j < box_rad ; ++j)
			for (k = -box_rad ; k < box_rad ; ++k) {
				val = model[(x+i)*size*size + (y+j)*size + (z+k)] ;
				avg += val ;
				rms += val*val ;
			}
			
			local_variation[x*size*size + y*size + z] = rms * weight - powf(avg * weight, 2.f) ;
			voxel_pos[(num_p++)*np + rank] = x*size*size + y*size + z ;
		}
		
		#pragma omp critical(num_vox)
		{
			if (num_p > num_vox)
				num_vox = num_p ;
		}
		if (rank == 0)
			num_vox *= np ;
	}
	
	// Sort and set keep highest num_supp values
	qsort(voxel_pos, num_vox, sizeof(long), compare_indices_variation) ;
	
	memset(supp, 0, vol*sizeof(uint8_t)) ;
	for (i = 0 ; i < num_supp ; ++i)
		supp[voxel_pos[i]] = 1 ;
}

void gen_rot(float rot[3][3], double *q) {
	double q0, q1, q2, q3, q01, q02, q03, q11, q12, q13, q22, q23, q33 ;
	
	q0 = q[0] ;
	q1 = q[1] ;
	q2 = q[2] ;
	q3 = q[3] ;
	
	q01 = q0*q1 ;
	q02 = q0*q2 ;
	q03 = q0*q3 ;
	q11 = q1*q1 ;
	q12 = q1*q2 ;
	q13 = q1*q3 ;
	q22 = q2*q2 ;
	q23 = q2*q3 ;
	q33 = q3*q3 ;
	
	rot[0][0] = (1. - 2.*(q22 + q33)) ;
	rot[0][1] = 2.*(q12 + q03) ;
	rot[0][2] = 2.*(q13 - q02) ;
	rot[1][0] = 2.*(q12 - q03) ;
	rot[1][1] = (1. - 2.*(q11 + q33)) ;
	rot[1][2] = 2.*(q01 + q23) ;
	rot[2][0] = 2.*(q02 + q13) ;
	rot[2][1] = 2.*(q23 - q01) ;
	rot[2][2] = (1. - 2.*(q11 + q22)) ;
}

/* Rotate and average intensity distribution using given quaternions and weights
 * Note: q=0 is at (0,0,0) and not (c,c,c)
 * Can set out = in
 */
void blur_intens(float *in, float *out) {
	#pragma omp parallel default(shared)
	{
		long r, i, j, x, y, z, vox[3], c = size / 2 ;
		float fx, fy, fz, cx, cy, cz, w ;
		float rot_vox[3], rot[3][3] ;
		float *priv_out = calloc(vol, sizeof(float)) ;
		
		#pragma omp for schedule(dynamic,1)
		for (r = 0 ; r < num_rot ; ++r) {
			if (quat[r*5 + 4] < 0.)
				continue ;
			
			gen_rot(rot, &(quat[r*5])) ;
			
			for (vox[0] = -c ; vox[0] < size-c ; ++vox[0])
			for (vox[1] = -c ; vox[1] < size-c ; ++vox[1])
			for (vox[2] = -c ; vox[2] < size-c ; ++vox[2]) {
				for (i = 0 ; i < 3 ; ++i) {
					rot_vox[i] = 0. ;
					
					for (j = 0 ; j < 3 ; ++j)
						rot_vox[i] += rot[i][j] * vox[j] ;
					
					rot_vox[i] = fmod(rot_vox[i] + size, size) ;
				}
				
				x = rot_vox[0] ;
				y = rot_vox[1] ;
				z = rot_vox[2] ;
				fx = rot_vox[0] - x ;
				fy = rot_vox[1] - y ;
				fz = rot_vox[2] - z ;
				cx = 1. - fx ;
				cy = 1. - fy ;
				cz = 1. - fz ;
				
				w = in[((vox[0]+size)%size)*size*size + ((vox[1]+size)%size)*size + ((vox[2]+size)%size)] * 
				    quat[r*5 + 4] ;
				
				priv_out[x*size*size + y*size + z]                                  += cx*cy*cz*w ;
				priv_out[x*size*size + y*size + ((z+1)%size)]                       += cx*cy*fz*w ;
				priv_out[x*size*size + ((y+1)%size)*size + z]                       += cx*fy*cz*w ;
				priv_out[x*size*size + ((y+1)%size)*size + ((z+1)%size)]            += cx*fy*fz*w ;
				priv_out[((x+1)%size)*size*size + y*size + z]                       += fx*cy*cz*w ;
				priv_out[((x+1)%size)*size*size + y*size + ((z+1)%size)]            += fx*cy*fz*w ;
				priv_out[((x+1)%size)*size*size + ((y+1)%size)*size + z]            += fx*fy*cz*w ;
				priv_out[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += fx*fy*fz*w ;
			}
		}
		
		#pragma omp barrier
		if (omp_get_thread_num() == 0)
			memset(out, 0, vol*sizeof(float)) ;
		#pragma omp barrier
		
		#pragma omp critical(priv_out)
		{
			for (x = 0 ; x < vol ; ++x)
				out[x] += priv_out[x] ;
		}
		
		free(priv_out) ;
	}
}

/* Match Bragg Fourier components
 * sigma parameter allows for small distance from input value at each voxel
 */
void match_bragg(fftwf_complex *fdens, float sigma) {
	long i ;
	fftwf_complex temp ;
	float mag ;
	
	if (sigma == 0.f) {
		for (i = 0 ; i < vol ; ++i)
		if (bragg_calc[i] != FLT_MAX)
			fdens[i] = bragg_calc[i] ;
	}
	else {
		for (i = 0 ; i < vol ; ++i) 
		if (bragg_calc[i] != FLT_MAX) {
			temp = fdens[i] - bragg_calc[i] ;
			mag = cabsf(temp) ;
			if (mag < sigma)
				fdens[i] = bragg_calc[i] ;
			else
				fdens[i] = bragg_calc[i] + sigma / mag * temp ;
		}
	}
}
		
