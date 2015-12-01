#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, char *argv[]) {
	long size, vol ;
	float *model, *pmodel ;
	double rescale ;
	FILE *fp ;
	
	if (argc < 5) {
		fprintf(stderr, "Format: %s <model_fname> <size> <rescale> <out_fname>\n", argv[0]) ;
		return 1 ;
	}
	size = atoi(argv[2]) ;
	rescale = atof(argv[3]) ;
	
	vol = size*size*size ;
	gsl_rng_env_setup() ;
	omp_set_num_threads(32) ;
	
	// Read model
	model = malloc(vol * sizeof(float)) ;
	fp = fopen(argv[1], "rb") ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Poisson sample rescale model
	pmodel = malloc(vol * sizeof(float)) ;
	#pragma omp parallel default(shared)
	{
		int rank = omp_get_thread_num() ;
		long x, y, z, c = size / 2 ;
		double val, dist, multip ;
		const gsl_rng_type *T ;
		gsl_rng *rng ;
		T = gsl_rng_default ;
		rng = gsl_rng_alloc(T) ;
		
		#pragma omp for schedule(static,1)
		for (x = 0 ; x < size ; ++x) {
			for (y = 0 ; y < size ; ++y)
			for (z = 0 ; z < size ; ++z) {
				dist = sqrt((x-c)*(x-c) + (y-c)*(y-c) + (z-c)*(z-c)) ;
				
				if (dist == 0.)
					multip = 293.33 ;
				else
					multip = 293.33 / dist ;
				
				val = model[x*size*size + y*size + z] * rescale * multip ;
				if (val > 100.)
					pmodel[x*size*size + y*size + z] = (val + gsl_ran_gaussian(rng, sqrt(val))) / multip ;
				else
					pmodel[x*size*size + y*size + z] = gsl_ran_poisson(rng, val) / multip ;
			}
			
			if (rank == 0)
				fprintf(stderr, "\r%ld/%ld", x, size) ;
		}
		
		gsl_rng_free(rng) ;
	}
	fprintf(stderr, "\n") ;
	
	// Save sampled model
	fp = fopen(argv[4], "wb") ;
	fwrite(pmodel, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	free(pmodel) ;
	
	return 0 ;
}
