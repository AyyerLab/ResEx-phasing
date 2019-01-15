#include "../../src/utils.h"
#include "../../src/fft.h"

long factorial(long i) {
	if (i < 2)
		return i ;
	if (i > 10) // Using Sterling's approximation anyway
		return -1 ;
	return i*factorial(i - 1) ;
}

int main(int argc, char *argv[]) {
	long x, y, z, size, c, vol, fac ;
	float vox_size, qscale, sigma, gamma ;
	double error, mean ;
	float *radius, *intens, *autocorr ;
	FILE *fp ;
	char fname[1024] ;
	struct fft_data fft ;
	
	if (argc < 5) {
		fprintf(stderr, "Liquidize: Apply liquid-like motion blurring to Fourier amplitudes\n") ;
		fprintf(stderr, "------------------------------------------------------------------\n") ;
		fprintf(stderr, "Needs three parameters:\n") ;
		fprintf(stderr, "\tres_at_edge: Resolution at edge of volume\n") ;
		fprintf(stderr, "\tsigma: RMS displacement of atom\n") ;
		fprintf(stderr, "\tgamma: Correlation length of atomic displacements\n") ;
		fprintf(stderr, "\nUsage: %s <cpx_fname> <res_at_edge> <sigma> <gamma>\n", argv[0]) ;
		fprintf(stderr, "<res_at_edge>, <sigma> and <gamma> must have consistent length units\n") ;
		fprintf(stderr, "\nOutput: <cpx_fname>-liq.raw\n") ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float complex)) ;
	c = size / 2 ;
	vol = size*size*size ;
	vox_size = atof(argv[2]) / 2. ;
	qscale = 1. / (2. * vox_size * c) ;
	sigma = atof(argv[3]) ;
	gamma = atof(argv[4]) ;

	// Allocate memory and FFTW plans
	fft_init(&fft, size, omp_get_max_threads()) ;
	fprintf(stderr, "Initialized fft\n") ;
	fft_create_plans(&fft) ;
	fprintf(stderr, "Created plans\n") ;
	autocorr = malloc(vol * sizeof(float)) ;
	intens = calloc(vol, sizeof(float)) ;
	radius = malloc(vol * sizeof(float)) ;
	
	// Read complex Fourier amplitudes
	fp = fopen(argv[1], "rb") ;
	fread(fft.rdensity, sizeof(float complex), vol, fp) ;
	fclose(fp) ;

	// Calculate intensities and shift origin
	mean = 0. ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
		  = powf(cabsf(fft.rdensity[x*size*size + y*size + z]), 2.f) ;
		radius[x*size*size + y*size + z] = sqrtf((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		mean += fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)] ;
	}

	// Get ideal autocorrelation
	fft_inverse(&fft) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z)
		autocorr[((x+c)%size)*size*size + ((y+c)%size)*size + ((z+c)%size)] = crealf(fft.rdensity[x*size*size + y*size + z]) / vol ;

	// Calculate modified intensities
	fprintf(stderr, "Starting convolutions\n") ;
	float *diff_array = malloc(vol * sizeof(float)) ;
	int break_flag = 0 ;
	#pragma omp parallel default(shared)
	{
		int n, omp_rank = omp_get_thread_num() ;
		long x, y, z, vox ;
		double diff, exponent ;
		for (n = 1 ; n <= 100 ; ++n) {
			if (break_flag)
				continue ;
			
			#pragma omp for schedule(static,1)
			for (x = 0 ; x < size ; ++x)
			for (y = 0 ; y < size ; ++y)
			for (z = 0 ; z < size ; ++z) {
				vox = x*size*size + y*size + z ;
				exponent = n * radius[vox] * vox_size / gamma ;
				fft.rdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]
				  = autocorr[vox] * expf(-exponent) ;
			}
			
			if (omp_rank == 0) {
				fft_forward(&fft) ;
				fac = factorial(n) ;
				error = 0. ;
			}
			#pragma omp barrier
			
			#pragma omp for schedule(static,1) reduction(+:error)
			for (x = 0 ; x < size ; ++x)
			for (y = 0 ; y < size ; ++y)
			for (z = 0 ; z < size ; ++z) {
				vox = x*size*size + y*size + z ;
				exponent = powf(2. * M_PI * sigma * radius[vox] * qscale, 2.f) ;
				if (n < 10) {
					diff = cabs(fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)]) *
						   pow(exponent, n) / fac * exp(-exponent) ;
				}
				else { // Stirling's approximation for n >= 10
					diff = log(cabs(fft.fdensity[((x+c+1)%size)*size*size + ((y+c+1)%size)*size + ((z+c+1)%size)])) +
						   n*log(exponent) - (log(2.*M_PI*n)/2 + n*log(n) - n) - exponent ;
					diff = exp(diff) ;
				}
				if (radius[vox] > 0.)
					diff /= 1. - exp(-exponent) ;
				diff_array[vox] = diff ;
				error += diff ;
				intens[vox] += diff ;
			}
			
			if (omp_rank == 0) {
				sprintf(fname, "data/liq-%.2d.raw", n) ;
				fp = fopen(fname, "wb") ;
				fwrite(diff_array, sizeof(float), vol, fp) ;
				fclose(fp) ;
				
				fprintf(stderr, "Finished n = %.3d, rel_error = %e\n", n, error/mean) ;
				if (n > 2 && error/mean < 1.e-5)
					break_flag = 1 ;
			}
			#pragma omp barrier
		}
	}

	// Save output
	sprintf(fname, "%s-liq.raw", remove_ext(argv[1])) ;
	fprintf(stderr, "Saving output to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;

	// Free memory
	fft_free(&fft) ;
	free(autocorr) ;
	free(intens) ;
	free(radius) ;
	free(diff_array) ;

	return 0 ;
}
