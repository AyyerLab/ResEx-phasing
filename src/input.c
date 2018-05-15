#include "input.h"

void input_init(struct input_data *self, long size) {
	self->size = size ;
	self->num_supp = 0 ;
	
	self->obs_mag = NULL ;
	self->bragg_calc = NULL ;
	self->support = NULL ;
	self->supp_locval = NULL ;
	self->inverse_cdf = NULL ;
	self->local_variation = NULL ;
}

int input_parse_intens(struct input_data *self, char *fname, float scale) {
	long i, j, k, s = self->size, c = s/2, vol = s*s*s ;
	float *intens ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found. Exiting.\n", fname) ;
		return 1 ;
	}
	
	intens = malloc(vol * sizeof(float)) ;
	self->obs_mag = malloc(vol * sizeof(float)) ;
	fread(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	fprintf(stderr, "Scale factor = %f\n", scale) ;
	
	for (i = 0 ; i < s ; ++i)
	for (j = 0 ; j < s ; ++j)
	for (k = 0 ; k < s ; ++k)
		if (intens[i*s*s + j*s + k] > 0.)
			self->obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= sqrt(intens[i*s*s + j*s + k]) * scale ;
		else
			self->obs_mag[((i+c+1)%s)*s*s + ((j+c+1)%s)*s + ((k+c+1)%s)]
				= intens[i*s*s + j*s + k] ;
	
	free(intens) ;
	
	return 0 ;
}

int input_parse_bragg(struct input_data *self, char *fname, float braggqmax) {
	long x, y, z, s = self->size, c = s/2, vol = s*s*s ;
	double distsq, c2 = c*c ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found. Exiting.\n", fname) ;
		return 1 ;
	}
	
	float complex *bragg_temp = malloc(vol * sizeof(float complex)) ;
	self->bragg_calc = malloc(vol * sizeof(float complex)) ;
	
	fread(bragg_temp, sizeof(float complex), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < s ; ++x)
	for (y = 0 ; y < s ; ++y)
	for (z = 0 ; z < s ; ++z) {
		distsq = ((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) / c2 ;
		
		if (distsq < braggqmax*braggqmax) 
			// Move (q=0) from center to corner
			self->bragg_calc[((x+c+1)%s)*s*s + ((y+c+1)%s)*s + ((z+c+1)%s)] 
			 = bragg_temp[x*s*s + y*s + z] 
			   * cexpf(-I * 2. * M_PI * (x-c+y-c+z-c) * c / s) ;
		else
			self->bragg_calc[((x+c+1)%s)*s*s + ((y+c+1)%s)*s + ((z+c+1)%s)] = FLT_MAX ;
	}
	
	free(bragg_temp) ;
	
	return 0 ;
}

int input_parse_support(struct input_data *self, char *fname) {
	long x, y, z, size = self->size, vol = size*size*size ;
	
	FILE *fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "%s not found.\n", fname) ;
		return 1 ;
	}
	self->support = malloc(vol * sizeof(uint8_t)) ;
	self->supp_locval = malloc(vol / 8 * sizeof(struct locval_pair)) ;
	fread(self->support, sizeof(uint8_t), vol, fp) ;
	fclose(fp) ;
	
	for (x = 0 ; x < 6 ; ++x)
		self->support_bounds[x] = size * (1 - (x % 2)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
	if (self->support[x*size*size + y*size + z]) {
		self->supp_locval[self->num_supp++].loc = x*size*size + y*size + z ;
		if (x < self->support_bounds[0])
			self->support_bounds[0] = x ;
		if (x > self->support_bounds[1])
			self->support_bounds[1] = x ;
		if (y < self->support_bounds[2])
			self->support_bounds[2] = y ;
		if (y > self->support_bounds[3])
			self->support_bounds[3] = y ;
		if (z < self->support_bounds[4])
			self->support_bounds[4] = z ;
		if (z > self->support_bounds[5])
			self->support_bounds[5] = z ;
	}
	
	fprintf(stderr, "num_supp = %ld\nSupport bounds: (", self->num_supp) ;
	for (x = 0 ; x < 6 ; ++x)
		fprintf(stderr, "%ld ", self->support_bounds[x]) ;
	fprintf(stderr, "\b)\n") ;
	
	return 0 ;
}

void input_init_iterate(struct input_data *self, char *fname, char *bg_fname, float *model, int do_bg_fitting) {
	long i, size = self->size, vol = size*size*size ;
	float val = sqrtf(vol) ;
	struct timeval t1 ;
	const gsl_rng_type *T ;
	gsl_rng *r ;
	FILE *fp ;
	int do_random_model = 0, do_init_bg = 0 ;
	
	fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Random start") ;
		do_random_model = 1 ;
	}
	else {
		fprintf(stderr, "Starting from %s", fname) ;
		fread(model, sizeof(float), vol, fp) ;
		fclose(fp) ;
	}
	
	if (do_bg_fitting) {
		fp = fopen(bg_fname, "rb") ;
		if (fp == NULL) {
			fprintf(stderr, " and with uniform background\n") ;
			do_init_bg = 1 ;
		}
		else {
			fprintf(stderr, " with background from %s\n", fname) ;
			fread(&(model[vol]), sizeof(float), vol, fp) ;
			fclose(fp) ;
		}
	}
	else {
		fprintf(stderr, "\n") ;
	}
	
	gsl_rng_env_setup() ;
	T = gsl_rng_default ;
	gettimeofday(&t1, NULL) ;
	r = gsl_rng_alloc(T) ;
	gsl_rng_set(r, t1.tv_sec + t1.tv_usec) ;
	
	if (do_random_model)
		memset(model, 0, vol*sizeof(float)) ;
	if (do_init_bg)
		memset(&(model[vol]), 0, vol*sizeof(float)) ;
	
	if (do_random_model || do_init_bg) {
		for (i = 0 ; i < vol ; ++i)
		if (do_random_model && self->support[i])
			model[i] = gsl_rng_uniform(r) ;
		
		if (do_bg_fitting) {
			for (i = 0 ; i < vol ; ++i)
			if (do_init_bg && self->obs_mag[i] > 0.)
				model[vol+i] = val ;
		}
	}
	
	gsl_rng_free(r) ;
}

int input_read_histogram(struct input_data *self, char *fname) {
	long m, i, num_hist ;
	double frac, sum_hist ;
	
	FILE *fp = fopen(fname, "r") ;
	if (fp == NULL) {
		fprintf(stderr, "Cannot find histogram file %s\n", fname) ;
		return 1 ;
	}
	fscanf(fp, "%ld", &num_hist) ;
	double *hist = malloc(num_hist * sizeof(double)) ;
	float *val = malloc(num_hist * sizeof(float)) ;
	double *cdf = malloc((num_hist+1) * sizeof(double)) ;
	for (i = 0 ; i < num_hist ; ++i)
		fscanf(fp, "%f %lf", &val[i], &hist[i]) ;
	fclose(fp) ;
	
	cdf[0] = 0. ;
	cdf[1] = hist[0] ;
	sum_hist = hist[0] ;
	for (i = 2 ; i <= num_hist ; ++i) {
		sum_hist += hist[i-1] ;
		cdf[i] = hist[i-1] + cdf[i-1] ;
	}
	for (i = 0 ; i <= num_hist ; ++i)
		cdf[i] /= sum_hist ;
	
	self->inverse_cdf = malloc(self->num_supp * sizeof(float)) ;
	self->inverse_cdf[0] = val[0] ;
	
	i = 1 ;
	for (m = 1 ; m < self->num_supp ; ++m) {
		frac = (double) m / self->num_supp ;
		if (i == num_hist - 1)
			self->inverse_cdf[m] = val[i] ;
		else if (frac <= cdf[i])
			self->inverse_cdf[m] = val[i-1] + (frac-cdf[i-1]) * (val[i] - val[i-1]) / (cdf[i] - cdf[i-1]) ;
		else {
			i++ ;
			while (frac > cdf[i])
				i++ ;
			self->inverse_cdf[m] = val[i-1] ;
		}
	}
	
	/*
	fp = fopen("data/inverse_cdf.raw", "wb") ; // Dump for debugging
	fwrite(self->inverse_cdf, sizeof(float), self->num_supp, fp) ;
	fclose(fp) ;
	*/
	
	free(hist) ;
	free(val) ;
	free(cdf) ;
	
	return 0 ;
}

int compare_indices(const void *a, const void *b) {
	struct locval_pair ia = *(struct locval_pair*) a ;
	struct locval_pair ib = *(struct locval_pair*) b ;
	return ia.val < ib.val ? -1 : ia.val > ib.val ;
}

/* Histogram matching
 * Matches histogram within support volume using inverse_cdf array.
 * Can be done in-place
 */
void input_match_histogram(struct input_data *self, float *in, float *out) {
	long i, vol = self->size*self->size*self->size ; 
	
	for (i = 0 ; i < self->num_supp ; ++i)
		self->supp_locval[i].val = in[self->supp_locval[i].loc] ;
	
	qsort(self->supp_locval, self->num_supp, sizeof(long), compare_indices) ;
	
	memset(out, 0, vol*sizeof(float)) ;
	for (i = 0 ; i < self->num_supp ; ++i)
		out[self->supp_locval[i].loc] = self->inverse_cdf[i] ;
}

/* Local variation support update
 * Similar to solvent flattening, keeping number of support voxels the same
 */
void input_update_support(struct input_data *self, float *model, long box_rad) {
	long i, num_vox = 0, size = self->size, vol = size*size*size ;
	
	if (self->local_variation == NULL)
		self->local_variation = malloc(vol * sizeof(struct locval_pair)) ;
	else
		memset(self->local_variation, 0, vol*sizeof(float)) ;
	#pragma omp parallel default(shared)
	{
		int np = omp_get_num_threads(), rank = omp_get_thread_num() ;
		long x, y, z, i, j, k, num_p = 0 ;
		float val, avg, rms, weight = 1. / powf(2 * box_rad + 1, 3.f) ;
		
		#pragma omp for schedule(static,1)
		for (x = self->support_bounds[0]-3*box_rad ; x <= self->support_bounds[1]+3*box_rad ; ++x)
		for (y = self->support_bounds[2]-3*box_rad ; y <= self->support_bounds[3]+3*box_rad ; ++y)
		for (z = self->support_bounds[4]-3*box_rad ; z <= self->support_bounds[5]+3*box_rad ; ++z) {
			avg = 0. ;
			rms = 0. ;
			for (i = -box_rad ; i < box_rad ; ++i)
			for (j = -box_rad ; j < box_rad ; ++j)
			for (k = -box_rad ; k < box_rad ; ++k) {
				val = model[(x+i)*size*size + (y+j)*size + (z+k)] ;
				avg += val ;
				rms += val*val ;
			}
			
			self->local_variation[num_p*np + rank].loc = x*size*size + y*size + z ;
			self->local_variation[num_p*np + rank].val = rms * weight - powf(avg * weight, 2.f) ;
			num_p++ ;
		}
		
		#pragma omp critical(num_vox)
		{
			if (num_p > num_vox)
				num_vox = num_p ;
		}
		#pragma omp barrier
		if (rank == 0)
			num_vox *= np ;
	}
	
	// Sort and set keep highest num_supp values
	qsort(self->local_variation, num_vox, sizeof(long), compare_indices) ;
	
	memset(self->support, 0, vol*sizeof(uint8_t)) ;
	for (i = 0 ; i < self->num_supp ; ++i)
		self->support[self->local_variation[i].loc] = 1 ;
}

/* Match Bragg Fourier components
 * sigma parameter allows for small distance from input value at each voxel
 */
void match_bragg(struct input_data *self, float complex *fdens, float delta) {
	long i, size = self->size, vol = size*size*size ;
	float complex temp ;
	float mag ;
	
	if (delta == 0.f) {
		for (i = 0 ; i < vol ; ++i)
		if (self->bragg_calc[i] != FLT_MAX)
			fdens[i] = self->bragg_calc[i] ;
	}
	else {
		for (i = 0 ; i < vol ; ++i) 
		if (self->bragg_calc[i] != FLT_MAX) {
			temp = fdens[i] - self->bragg_calc[i] ;
			mag = cabsf(temp) ;
			if (mag < delta)
				fdens[i] = self->bragg_calc[i] ;
			else
				fdens[i] = self->bragg_calc[i] + delta / mag * temp ;
		}
	}
}

void input_free(struct input_data *self) {
	if (self->obs_mag != NULL)
		free(self->obs_mag) ;
	if (self->bragg_calc != NULL)
		free(self->bragg_calc) ;
	if (self->support != NULL)
		free(self->support) ;
	
	// Histogram matching
	if (self->supp_locval != NULL)
		free(self->supp_locval) ;
	if (self->inverse_cdf != NULL)
		free(self->inverse_cdf) ;
	
	// Local variation support update
	if (self->local_variation != NULL)
		free(self->local_variation) ;
}
