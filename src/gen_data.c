#include "brcont.h"

float gaussian(float x, float width) {
	return exp(- x*x / 2 / width/width) ;
}

int main(int argc, char *argv[]) {
	long x, y, z, center, num_bins, num_cells ;
	long spotnum[3], i[3], ksize, krad ;
	float dist, sigma, blur ;
	float *dwfactor, *kernel ;
	long a[3] = {3, 4, 7} ;
	float q_scale[3] = {1.1639, 1.09405, 1.1691} ;
	
	long t[3] ;
	FILE *fp ;
	
	omp_set_num_threads(32) ;
	sigma = 1.5 ;
	blur = 0.2 ;
	num_cells = 1000 ;
	
	// Setup arrays and parse files
	if (setup_gen())
		return 2 ;
	fprintf(stderr, "Finished setup\n") ;
	
	center = size / 2 ;
	
	// Generate q-dependent factor
	num_bins = (int) size ;
	dwfactor = malloc(num_bins * sizeof(float)) ;
	for (x = 0 ; x < num_bins ; ++x)
		dwfactor[x] = expf(-powf(2.*M_PI*1.e-3*sigma*x, 2.)) ;
	
	// Generate lattice function
	ksize = (int) 10 * blur ;
	krad = (ksize - 1) / 2 ;
	kernel = malloc(ksize * ksize * ksize * sizeof(float)) ;
	for (i[0] = 0 ; i[0] < ksize ; ++i[0])
	for (i[1] = 0 ; i[1] < ksize ; ++i[1])
	for (i[2] = 0 ; i[2] < ksize ; ++i[2])
		kernel[i[0]*ksize*ksize + i[1]*ksize + i[2]] = 
			gaussian(
				sqrt((krad-i[0])*(krad-i[0]) + 
					(krad-i[1])*(krad-i[1]) + 
					(krad-i[2])*(krad-i[2])), 
				blur) ;
	fprintf(stderr, "Generated kernel with ksize = %ld\n", ksize) ;
	
	for (x = 0 ; x < 3 ; ++x)
		spotnum[x] = (int) floor((size-ksize) / a[x] / 2.) ;
	fprintf(stderr, "spotnums = %ld, %ld, %ld\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
	for (i[0] = -spotnum[0] ; i[0] <= spotnum[0] ; ++i[0])
	for (i[1] = -spotnum[1] ; i[1] <= spotnum[1] ; ++i[1])
	for (i[2] = -spotnum[2] ; i[2] <= spotnum[2] ; ++i[2]) {
		for (x = 0 ; x < 3 ; ++x)
			t[x] = a[x] * i[x] + center ;
		
		// Convolve with Gaussian spot kernel
		for (x = 0 ; x < ksize ; ++x)
		for (y = 0 ; y < ksize ; ++y)
		for (z = 0 ; z < ksize ; ++z)
			intens[(t[0]+x-krad)*size*size + (t[1]+y-krad)*size + t[2]+z-krad] 
			 = kernel[x*ksize*ksize + y*ksize + z] ;
	}
	fprintf(stderr, "Generated object with spotnums = %ld, %ld, %ld\n", spotnum[0], spotnum[1], spotnum[2]) ;
	
	// Symmetrize intensity
	symmetrize_incoherent(bragg_calc, incoh_mag) ;
	symmetrize_coherent(bragg_calc, coh_mag) ;
	fprintf(stderr, "Symmetrized intensity\n") ;
	
	// Apply q-dependent modulation
	#pragma omp parallel for default(shared) private(x,y,z,dist) schedule(static,1)
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dist = pow((x-center)*q_scale[0], 2.) ;
		dist += pow((y-center)*q_scale[1], 2.) ;
		dist += pow((z-center)*q_scale[2], 2.) ;
		dist = sqrt(dist) ;
		
		intens[x*size*size + y*size + z]
		 *= num_cells * powf(/*dwfactor[(int) dist] **/ 
//		        coh_mag[((x+center)%size)*size*size + ((y+center)%size)*size + ((z+center)%size)], 2.f) ;
		        coh_mag[x*size*size + y*size + z], 2.f) ;
		
		intens[x*size*size + y*size + z] 
//		 += powf(/*(1. - dwfactor[(int) dist]) **/
		 = powf(/*(1. - dwfactor[(int) dist]) **/
//		        incoh_mag[x*size*size + y*size + z], 2.f) ;
		        cabsf(bragg_calc[((x+center)%size)*size*size + ((y+center)%size)*size + ((z+center)%size)]), 2.f) ;
	}
	fprintf(stderr, "Applied DW factor\n") ;
	
	// Write to file
	fp = fopen("data/3wu2_dimer-intens.raw", "wb") ;
	fwrite(intens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	return 0 ;
}
