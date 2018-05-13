#include "fft.h"

void fft_init(struct fft_data *self, long size, int num_threads) {
	self->size = size ;
	long vol = size*size*size ;
	
	fftwf_init_threads() ;
	fftwf_plan_with_nthreads(num_threads) ;
	sprintf(self->wisdom_fname, "data/wisdom_%ld_%d", size, num_threads) ;
	
	self->rdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
	self->fdensity = fftwf_malloc(vol * sizeof(fftwf_complex)) ;
}

void fft_create_plans(struct fft_data *self) {
	long size = self->size ;
	FILE *fp = fopen(self->wisdom_fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Measuring plans\n") ;
		self->forward_plan = fftwf_plan_dft_3d(size, size, size, self->rdensity, self->fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		self->inverse_plan = fftwf_plan_dft_3d(size, size, size, self->fdensity, self->rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
		
		fp = fopen(self->wisdom_fname, "wb") ;
		fftwf_export_wisdom_to_file(fp) ;
		fclose(fp) ;
	
		fprintf(stderr, "Created plans. Saved to %s\n", self->wisdom_fname) ;
	}
	else {
		fftwf_import_wisdom_from_file(fp) ;
		fclose(fp) ;
		
		self->forward_plan = fftwf_plan_dft_3d(size, size, size, self->rdensity, self->fdensity, FFTW_FORWARD, FFTW_MEASURE) ;
		self->inverse_plan = fftwf_plan_dft_3d(size, size, size, self->fdensity, self->rdensity, FFTW_BACKWARD, FFTW_MEASURE) ;
	}
}

void fft_gaussian_blur(struct fft_data *self, float *model, float blur) {
	long x, y, z, size = self->size, c = size/2, vol = size*size*size ;
	float rsq, fblur = size / (2. * M_PI * blur) ;
	
	for (x = 0 ; x < vol ; ++x)
		self->rdensity[x] = model[x] ;
	
	// Blur density
	fftwf_execute(self->forward_plan) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		rsq = ((x+c)%size-c)*((x+c)%size-c) + ((y+c)%size-c)*((y+c)%size-c) + ((z+c)%size-c)*((z+c)%size-c) ;
		self->fdensity[x*size*size + y*size + z] *= expf(-rsq / 2. / fblur / fblur) ;
	}
	
	fftwf_execute(self->inverse_plan) ;
	
	for (x = 0 ; x < vol ; ++x)
		model[x] = crealf(self->rdensity[x]) / vol ;
}

void fft_apply_shrinkwrap(struct fft_data *self, float *model, float blur, float threshold, uint8_t *support, char *fname) {
	long x, size = self->size, vol = size*size*size ;
	FILE *fp ;
	
	fft_gaussian_blur(self, model, blur) ;
	
	// Apply threshold
	for (x = 0 ; x < vol ; ++x)
	if (model[x] > threshold)
		support[x] = 1 ;
	else
		support[x] = 0 ;
	
	if (fname != NULL) {
		mkdir(dirname(fname), S_IRWXU|S_IRGRP|S_IROTH) ;
		fp = fopen(fname, "wb") ;
		fwrite(support, sizeof(uint8_t), vol, fp) ;
		fclose(fp) ;
		
		//fp = fopen("data/smoothed.raw", "wb") ;
		//fwrite(model, sizeof(float), vol, fp) ;
		//fclose(fp) ;
	}
}

void fft_free(struct fft_data* self) {
	fftwf_destroy_plan(self->forward_plan) ;
	fftwf_destroy_plan(self->inverse_plan) ;
}
