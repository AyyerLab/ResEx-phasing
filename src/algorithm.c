#include "algorithm.h"

/* Fourier-space projection 
 * Makes model consistent with data. Modifies Fourier magnitudes and keeps
 * phases unchanged. Applies symmetry depending on the specified point group.
 * Can be applied in-place (out = in)
 */
static void proj_fourier(struct algorithm_data *self, float *in, float *out) {
	long i, vol = self->vol ;
	float norm_factor = 1. / (float) self->vol ;
	struct volume_data *volume = self->volume ;
	struct input_data *input = self->input ;
	struct fft_data *fft = self->fft ;
	
	if (self->do_bg_fitting) {
		for (i = 0 ; i < vol ; ++i) {
			fft->rdensity[i] = in[i] ;
			out[vol+i] = in[vol+i] ;
		}
	}
	else {
		for (i = 0 ; i < vol ; ++i)
			fft->rdensity[i] = in[i] ;
	}
	
	fft_forward(fft) ;
	
	if (self->do_bg_fitting)
		volume_symmetrize_incoherent(volume, fft->fdensity, self->exp_mag, &(out[vol])) ;
	else
		volume_symmetrize_incoherent(volume, fft->fdensity, self->exp_mag, NULL) ;
	
	if (self->do_blurring)
		volume_rotational_blur(volume, self->exp_mag, self->exp_mag, self->quat) ;
	
	match_bragg(input, fft->fdensity, 0.f) ;
	
	if (self->do_bg_fitting) {
		for (i = 0 ; i < vol ; ++i) {
			if (input->obs_mag[i] > 0.) {
				fft->fdensity[i] *= input->obs_mag[i] / self->exp_mag[i] ;
				out[vol+i] = in[vol+i] * input->obs_mag[i] / self->exp_mag[i] ;
			}
			else if (input->obs_mag[i] == 0.) {
				fft->fdensity[i] = 0. ;
				out[vol+i] = 0. ;
			}
			else {
				out[vol+i] = 0. ;
			}
		}
	}
	else {
		for (i = 0 ; i < vol ; ++i) {
			if (input->obs_mag[i] > 0.)
				fft->fdensity[i] *= input->obs_mag[i] / self->exp_mag[i] ;
			else if (input->obs_mag[i] == 0.)
				fft->fdensity[i] = 0. ;
		}
	}
	
	fft_inverse(fft) ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = crealf(fft->rdensity[i]) * norm_factor ;
}

/* Direct-space projection
 * Applies the finite support constraint in direct/real-space
 * In addition one can apply the following additional constraints by
 * setting the following global variables.
 * 	do_bg_fitting - Azimuthally average background part of iterate
 * 	do_positivity - Set negative values inside support to zero
 * 	do_histogram - Project values inside support such that they match a 
 * 	               target histogram.
 * 	do_local_variation - Calculate local variation and update support keeping
 * 	                     support size the same
 * Can be applied in-place (out = in)
 */
static void proj_direct(struct algorithm_data *self, float *in, float *out) {
	long vol = self->vol ;
	struct volume_data *volume = self->volume ;
	struct input_data *input = self->input ;
	
	if (self->do_local_variation > 0)
		input_update_support(input, in, 2) ;
	
	if (self->do_histogram) {
		input_match_histogram(input, in, out) ;
	}
	else {
		long i ;
		for (i = 0 ; i < vol ; ++i)
			out[i] = in[i] * input->support[i] ;
		
		if (self->do_positivity) {
			for (i = 0 ; i < vol ; ++i)
			if (out[i] < 0.)
				out[i] = 0. ;
		}
	}
	
	if (self->do_bg_fitting)
		volume_radial_average(volume, &(in[vol]), &(out[vol])) ;
}

/* Difference Map algorithm
 * Update rule (LaTeX syntax)
 * x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right)P_F(x_n) - \frac{x_n}{\beta}\right] - P_F\left[\left(1-\frac{1}{\beta}\right)P_D(x_n) + \frac{x_n}{\beta}\right]\right\}
 * Same as all other algorithms for beta = 1
 */
float DM_algorithm(struct algorithm_data *self) {
	long i ;
	float change = 0. ;
	
	proj_fourier(self, self->iterate, self->p1) ;
	if (self->beta != 1.)
		proj_direct(self, self->iterate, self->p2) ;
	
	for (i = 0 ; i < self->num_vox ; ++i) {
		self->r1[i] = (1. + 1./self->beta) * self->p1[i] - self->iterate[i] / self->beta ;
		if (self->beta != 1.)
			self->r2[i] = (1. - 1./self->beta) * self->p2[i] + self->iterate[i] / self->beta ;
	}
	
	proj_direct(self, self->r1, self->p2) ;
	if (self->beta != 1.)
		proj_fourier(self, self->r2, self->p1) ;
	
	for (i = 0 ; i < self->num_vox ; ++i) {
		float diff = self->beta * (self->p2[i] - self->p1[i]) ;
		self->iterate[i] += diff ;
		if (i < self->vol)
			change += diff*diff ;
	}
	
	return sqrtf(change / self->vol) ;
}

/* Modified Difference Map algorithm
 * Update rule (LaTeX syntax)
 * x_n' = \beta x_n + (1-\beta) P_F(x_n)
 * x_{n+1} = x_n' + P_F\left[2 P_D(x_n') - x_n'\right] - P_D(x_n')
 * Same as all other algorithms for beta = 1
 */
float mod_DM_algorithm(struct algorithm_data *self) {
	long i ;
	float change = 0. ;
	
	if (self->beta != 1.) {
		proj_fourier(self, self->iterate, self->p1) ;
		
		for (i = 0 ; i < self->num_vox ; ++i)
			self->iterate[i] = self->beta*self->iterate[i] + (1. - self->beta) * self->p1[i] ;
	}
	
	proj_direct(self, self->iterate, self->p2) ;
	
	for (i = 0 ; i < self->num_vox ; ++i)
		self->r1[i] = 2. * self->p2[i] - self->iterate[i] ;
	
	proj_fourier(self, self->r1, self->p1) ;
	
	for (i = 0 ; i < self->num_vox ; ++i) {
		float diff = self->p1[i] - self->p2[i] ;
		self->iterate[i] += diff ;
		if (i < self->vol)
			change += diff*diff ;
	}
	
	return sqrtf(change / self->vol) ;
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
float RAAR_algorithm(struct algorithm_data *self) {
	long i ;
	float change = 0. ;
	
	proj_fourier(self, self->iterate, self->p1) ;
	
	for (i = 0 ; i < self->num_vox ; ++i)
		self->r1[i] = 2. * self->p1[i] ;
	proj_direct(self, self->r1, self->r2) ;
	
	for (i = 0 ; i < self->num_vox ; ++i)
		self->r1[i] = - self->p1[i] ;
	proj_direct(self, self->r1, self->p2) ;
	
	for (i = 0 ; i < self->num_vox ; ++i) {
		float diff = (self->beta - 1.) * self->iterate[i] + self->beta * (self->r2[i] + self->p2[i]) + (1. - 2. * self->beta) * self->p1[i] ;
		self->iterate[i] += diff ;
		if (i < self->vol)
			change += diff*diff ;
	}
	
	return sqrtf(change / self->vol) ;
}

/* HIO algorithm
 * Update rule (LaTeX syntax)
 * x_{n+1} = x_n + \beta \left\{P_D\left[\left(1+\frac{1}{\beta}\right) P_F(x_n) - \frac{x_n}{\beta}\right] - P_F(x_n)\right\}
 * Same as all other algorithms for beta = 1
 */
float HIO_algorithm(struct algorithm_data *self) {
	long i ;
	float change = 0. ;
	
	proj_fourier(self, self->iterate, self->p1) ;
	
	for (i = 0 ; i < self->num_vox ; ++i)
		self->r1[i] = (1. + 1./self->beta) * self->p1[i] - self->iterate[i] / self->beta ;
	
	proj_direct(self, self->r1, self->p2) ;
	
	for (i = 0 ; i < self->num_vox ; ++i) {
		float diff = self->beta * (self->p2[i] - self->p1[i]) ;
		self->iterate[i] += diff ;
		if (i < self->vol)
			change += diff*diff ;
	}
	
	return sqrtf(change / self->vol) ;
}

/* Error Reduction algorithm
 * Update rule (LaTeX style)
 * x_{n+1} = P_D[P_F(x_n)]
 * Obviously different from others. Use only in averaging phase.
 */
float ER_algorithm(struct algorithm_data *self) {
	long i ;
	float change = 0. ;
	
	proj_fourier(self, self->iterate, self->p1) ;
	proj_direct(self, self->p1, self->p2) ;
	
	for (i = 0 ; i < self->num_vox ; ++i) {
		self->iterate[i] = self->p2[i] ;
		if (i < self->vol) {
			float diff = self->p2[i] - self->p1[i] ;
			change += diff*diff ;
		}
	}
	
	return sqrtf(change / self->vol) ;
}

/* ============================================================ */

void make_recon_folders(struct algorithm_data *self) {
	char fname[1024] ;
	
	sprintf(fname, "%s-slices", self->output_prefix) ;
	mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	sprintf(fname, "%s-fslices", self->output_prefix) ;
	mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	if (self->do_bg_fitting) {
		sprintf(fname, "%s-radavg", self->output_prefix) ;
		mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	}
	if (self->do_local_variation) {
		sprintf(fname, "%s-support", self->output_prefix) ;
		mkdir(fname, S_IRWXU|S_IRGRP|S_IROTH) ;
	}
}

float run_iteration(struct algorithm_data *self, int iter) {
	char algo[8] ;
	float error ;
	
	if (iter <= self->num_iter)
		strcpy(algo, self->algorithms[iter-1]) ;
	else
		strcpy(algo, self->avg_algorithms[iter-self->num_iter-1]) ;
	
	if (strcmp(algo, "DM") == 0)
		error = DM_algorithm(self) ;
	else if (strcmp(algo, "HIO") == 0)
		error = HIO_algorithm(self) ;
	else if (strcmp(algo, "RAAR") == 0)
		error = RAAR_algorithm(self) ;
	else if (strcmp(algo, "mod-DM") == 0)
		error = mod_DM_algorithm(self) ;
	else if (strcmp(algo, "ER") == 0)
		error = ER_algorithm(self) ;
	else {
		fprintf(stderr, "Could not understand algorithm name: %s\n", algo) ;
		return -1 ;
	}
	
	return error ;
}

void save_current(struct algorithm_data *self, int iter, struct timeval t1, struct timeval t2, float error) {
	char fname[1024] ;
	FILE *fp ;
	
	sprintf(fname, "%s-log.dat", self->output_prefix) ;
	fp = fopen(fname, "a") ;
	fprintf(fp, "%.4d\t%.2f s\t%f\n", 
			iter, (double)(t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000., error) ;
	fclose(fp) ;
	
	sprintf(fname, "%s-slices/%.4d.raw", self->output_prefix, iter) ;
	volume_dump_slices(self->p1, fname, self->size, 0) ;
	sprintf(fname, "%s-fslices/%.4d.raw", self->output_prefix, iter) ;
	volume_dump_slices(self->exp_mag, fname, self->size, 1) ;
	if (self->do_local_variation) {
		sprintf(fname, "%s-support/%.4d.supp", self->output_prefix, iter) ;
		volume_dump_support_slices(self->input->support, fname, self->size) ;
	}
	if (self->do_bg_fitting) {
		sprintf(fname, "%s-radavg/%.4d.raw", self->output_prefix, iter) ;
		fp = fopen(fname, "wb") ;
		fwrite(self->volume->radavg, sizeof(double), self->size/2, fp) ;
		fclose(fp) ;
	}
	
	if (iter == self->num_iter) {
		sprintf(fname, "%s-last.raw", self->output_prefix) ;
		fp = fopen(fname, "w") ;
		fwrite(self->iterate, sizeof(float), self->num_vox, fp) ;
		fclose(fp) ;
	}
}

void save_output(struct algorithm_data *self) {
	char fname[1024] ;
	FILE *fp ;
	
	sprintf(fname, "%s-pf.raw", self->output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(self->average_p1, sizeof(float), self->vol, fp) ;
	fclose(fp) ;
	
	sprintf(fname, "%s-pd.raw", self->output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(self->average_p2, sizeof(float), self->vol, fp) ;
	fclose(fp) ;
	
	if (self->do_bg_fitting) {
		sprintf(fname, "%s-bg.raw", self->output_prefix) ;
		fp = fopen(fname, "wb") ;
		fwrite(&(self->p1[self->vol]), sizeof(float), self->vol, fp) ;
		fclose(fp) ;
		
		sprintf(fname, "%s-radavg.raw", self->output_prefix) ;
		fp = fopen(fname, "wb") ;
		fwrite(self->volume->radavg, sizeof(float), self->size/2, fp) ;
		fclose(fp) ;
	}
	
	if (self->do_local_variation) {
		sprintf(fname, "%s-supp.supp", self->output_prefix) ;
		fp = fopen(fname, "wb") ;
		fwrite(self->input->support, sizeof(uint8_t), self->vol, fp) ;
		fclose(fp) ;
	}
}

static int test_alg_name(char *name) {
	if (strcmp(name, "DM") != 0 &&
	    strcmp(name, "HIO") != 0 &&
	    strcmp(name, "RAAR") != 0 &&
	    strcmp(name, "mod-DM") != 0 &&
	    strcmp(name, "ER") != 0)
		return 1 ;
	else
		return 0 ;
}

/* Generates list of algorithms given algorithm and avg_algorithm strings.
 */
int parse_algorithm_strings(struct algorithm_data *self, char *alg_string, char *avg_string) {
	char *token, string[9999] ;
	int i, current_num = 0, is_num = 1 ;
	self->num_iter = 0 ;
	self->num_avg_iter = 0 ;
	
	strcpy(string, alg_string) ;
	token = strtok(string, " \n") ;
	while (token != NULL) {
		if (is_num) {
			self->num_iter += atoi(token) ;
			is_num = 0 ;
		}
		else {
			if (test_alg_name(token)) {
				fprintf(stderr, "Unknown algorithm %s in algorithm %s\n", token, alg_string) ;
				return 1 ;
			}
			is_num = 1 ;
		}
		token = strtok(NULL, " \n") ;
	}
	fprintf(stderr, "Total number of normal iterations = %d\n", self->num_iter) ;
	
	self->algorithms = malloc(self->num_iter * sizeof(*(self->algorithms))) ;
	is_num = 1 ;
	self->num_iter = 0 ;
	strcpy(string, alg_string) ;
	token = strtok(string, " \n") ;
	while (token != NULL) {
		if (is_num) {
			current_num = atoi(token) ;
			is_num = 0 ;
		}
		else {
			for (i = self->num_iter ; i < self->num_iter + current_num ; ++i)
				strcpy(self->algorithms[i], token) ;
			self->num_iter += current_num ;
			is_num = 1 ;
		}
		token = strtok(NULL, " \n") ;
	}
	
	strcpy(string, avg_string) ;
	token = strtok(string, " \n") ;
	while (token != NULL) {
		if (is_num) {
			self->num_avg_iter += atoi(token) ;
			is_num = 0 ;
		}
		else {
			if (test_alg_name(token)) {
				fprintf(stderr, "Unknown algorithm %s in avg_algorithm %s\n", token, avg_string) ;
				return 1 ;
			}
			is_num = 1 ;
		}
		token = strtok(NULL, " \n") ;
	}
	fprintf(stderr, "Total number of averaging iterations = %d\n", self->num_avg_iter) ;
	
	self->avg_algorithms = malloc(self->num_avg_iter * sizeof(*(self->avg_algorithms))) ;
	is_num = 1 ;
	self->num_avg_iter = 0 ;
	strcpy(string, avg_string) ;
	token = strtok(string, " \n") ;
	while (token != NULL) {
		if (is_num) {
			current_num = atoi(token) ;
			is_num = 0 ;
		}
		else {
			for (i = self->num_avg_iter ; i < self->num_avg_iter + current_num ; ++i)
				strcpy(self->avg_algorithms[i], token) ;
			self->num_avg_iter += current_num ;
			is_num = 1 ;
		}
		token = strtok(NULL, " \n") ;
	}
	
	return 0 ;
}

void algorithm_allocate_memory(struct algorithm_data *self) {
	self->iterate = malloc(self->num_vox * sizeof(float)) ;
	self->exp_mag = malloc(self->vol * sizeof(float)) ;
	
	self->p1 = calloc(self->num_vox, sizeof(float)) ;
	self->p2 = malloc(self->num_vox * sizeof(float)) ;
	self->average_p1 = calloc(self->num_vox, sizeof(float)) ;
	self->average_p2 = malloc(self->num_vox * sizeof(float)) ;
	self->r1 = malloc(self->num_vox * sizeof(float)) ;
	if (self->beta != 1.)
		self->r2 = malloc(self->num_vox * sizeof(float)) ;
}

void calc_prtf(struct algorithm_data *self, int num_bins) {
	long x, y, z, bin ;
	long dx, dy, dz, size = self->size, vol = self->vol, c = size/2, c1 = size/2+1 ;
	float obs_val, *prtf = calloc(num_bins, sizeof(float)) ;
	long *bin_count = calloc(num_bins, sizeof(long)) ;
	FILE *fp ;
	char fname[1024] ;
	float *model1 = self->average_p1 ;
	float *model2 = self->average_p2 ;
	struct volume_data *volume = self->volume ;
	struct input_data *input = self->input ;
	struct fft_data *fft = self->fft ;
	
	// FFT average p2 model
	for (x = 0 ; x < vol ; ++x)
		fft->rdensity[x] = model2[x] ;
	
	fft_forward(fft) ;
	
	// Calculate exp_mag for average p2 model
	if (self->do_bg_fitting)
		volume_symmetrize_incoherent(volume, fft->fdensity, self->exp_mag, &(model2[vol])) ;
	else
		volume_symmetrize_incoherent(volume, fft->fdensity, self->exp_mag, NULL) ;
	
	// Save exp_mag
	sprintf(fname, "%s-expmag.raw", self->output_prefix) ;
	volume_dump_slices(self->exp_mag, fname, size, 1) ;
	
	// Calculate PRTF by comparing with obs_mag
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = (x + c1) % size - c1 ;
		dy = (y + c1) % size - c1 ;
		dz = (z + c1) % size - c1 ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / c1 * num_bins + 0.5 ;
		obs_val = input->obs_mag[x*size*size + y*size + z] ;
		
		if (bin < num_bins && obs_val > 0.) {
			prtf[bin] += self->exp_mag[x*size*size + y*size + z] / obs_val ;
			bin_count[bin]++ ;
		}
		else if (bin < num_bins && obs_val == 0.)
			bin_count[bin]++ ;
		
		self->p2[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)]
			= pow(cabsf(fft->fdensity[x*size*size + y*size + z]), 2.) ;
	}
	
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (bin_count[bin] > 0)
		prtf[bin] /= bin_count[bin] ;
	else 
		prtf[bin] = 1. ;
	
	// Save PRTF
	sprintf(fname, "%s-prtf.dat", self->output_prefix) ;
	fp = fopen(fname, "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, prtf[bin]) ;
	fclose(fp) ;
	
	// Save frecon (intensity from model)
	sprintf(fname, "%s-frecon.raw", self->output_prefix) ;
	fp = fopen(fname, "wb") ;
	fwrite(self->p2, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// If needed, normalize models by PRTF
	if (self->do_normalize_prtf) {
		fprintf(stderr, "Normalizing p2 model by PRTF\n") ;
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			dx = (x + c1) % size - c1 ;
			dy = (y + c1) % size - c1 ;
			dz = (z + c1) % size - c1 ;
			bin = sqrt(dx*dx + dy*dy + dz*dz) / c1 * num_bins + 0.5 ;
			
			if (bin < num_bins && prtf[bin] > 0.)
				fft->fdensity[x*size*size + y*size + z] /= prtf[bin] ;
			else
				fft->fdensity[x*size*size + y*size + z] = 0. ;
		}
		
		fft_inverse(fft) ;
		
		for (x = 0 ; x < vol ; ++x) {
			model2[x] = crealf(fft->rdensity[x]) / vol ;
			fft->rdensity[x] = model1[x] ;
		}
		fft_forward(fft) ;
		
		fprintf(stderr, "Normalizing p1 model by PRTF\n") ;
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			dx = (x + c1) % size - c1 ;
			dy = (y + c1) % size - c1 ;
			dz = (z + c1) % size - c1 ;
			bin = sqrt(dx*dx + dy*dy + dz*dz) / c1 * num_bins + 0.5 ;
			
			if (bin < num_bins && prtf[bin] > 0.)
				fft->fdensity[x*size*size + y*size + z] /= prtf[bin] ;
			else
				fft->fdensity[x*size*size + y*size + z] = 0. ;
		}
		
		fft_inverse(fft) ;
		for (x = 0 ; x < vol ; ++x)
			model1[x] = crealf(fft->rdensity[x]) / vol ;
	}
	
	free(prtf) ;
	free(bin_count) ;
}

