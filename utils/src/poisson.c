#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../../src/utils.h"

int main(int argc, char *argv[]) {
	long i, x, y, z, size, c, vol, num_vox = 0 ;
	float *model, *pmodel ;
	double rescale, numr, denr, target, dist ;
	FILE *fp ;
	
	if (argc < 4) {
		fprintf(stderr, "Format: %s <model_fname> <target_count_at_edge> <out_fname>\n", argv[0]) ;
		return 1 ;
	}
	size = get_size(argv[1], sizeof(float)) ;
	target = atof(argv[2]) ;
	
	c = size/2 ;
	vol = size*size*size ;
	gsl_rng_env_setup() ;
	omp_set_num_threads(32) ;
	
	// Read model
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Rescale model to match target_count_at_edge
	rescale = 0. ;
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dist = sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
		if ((int) dist == c-1) {
			rescale += model[x*size*size + y*size + z] ;
			num_vox++ ;
		}
	}
	
	rescale = target / (rescale / num_vox) ;
	
	// Poisson sample rescaled model
	pmodel = malloc(vol * sizeof(float)) ;
	#pragma omp parallel default(shared)
	{
		int rank = omp_get_thread_num() ;
		long vox, x, y, z, c = size / 2 ;
		double val, dist, multip ;
		const gsl_rng_type *T ;
		gsl_rng *rng ;
		T = gsl_rng_default ;
		rng = gsl_rng_alloc(T) ;
		
		#pragma omp for schedule(static,1)
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			vox = x*size*size + y*size + z ;
			dist = sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
			
			if (dist == 0.)
				multip = c ;
			else
				multip = c / dist ;
			
			val = model[vox] * rescale * multip ;
			if (val > 100.)
				pmodel[vox] = (val + gsl_ran_gaussian(rng, sqrt(val))) / multip ;
			else if (val > 0.)
				pmodel[vox] = gsl_ran_poisson(rng, val) / multip ;
			else
				pmodel[vox] = model[vox] ;
			
			if (model[vox] > 0.) {
				numr += model[vox] * pmodel[vox] ;
				denr += pmodel[vox] * pmodel[vox] ;
			}
			
			if (rank == 0 && y == size - 1 && z == size - 1)
				fprintf(stderr, "\rFinished x = %ld", x) ;
		}
		
		gsl_rng_free(rng) ;
	}
	
	rescale = numr / denr ;
	fprintf(stderr, "\nScale factor: %e\n", rescale) ;
	for (i = 0 ; i < vol ; ++i)
		pmodel[i] *= rescale ;
	
	// Save sampled model
	fp = fopen(argv[3], "wb") ;
	fwrite(pmodel, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	free(pmodel) ;
	
	return 0 ;
}
