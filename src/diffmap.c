/* Phasing both Bragg and continuous intensities using a divide and 
	concur approach. Both contraints are applied to different models
	in the divide constraint. The third constraint is the support
	constraint which is applied in real-space. The concur constraint 
	enforces equality between the three models.
	
	The iterate is a list of three vectors, each of size 'vol' i.e.
		Dimensions[x] = {3, vol}
*/

#include "brcont.h"

// Symmetrize intensity according to P(2_1)(2_1)(2_1) space group
// (size,hsize,ksize,lsize)
void symmetrize_intens(fftwf_complex *in, float *out, int flag) {
	long hs, ks, ls, kc, lc ;
	
	if (flag == 1) {
		hs = size ;
		ks = size ;
		ls = size ;
		kc = ks / 2 ;
		lc = ls / 2 ;
	}
	else if (flag == 0) {
		hs = hsize ;
		ks = ksize ;
		ls = lsize ;
		kc = ks / 2 ;
		lc = ls / 2 ;
	}
	
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
			}
			
			for (z = 1 ; z <= lc ; ++z) {
				ave_intens = 0.5 * (
					powf(cabsf(in[x*ks*ls + z]), 2.) + 
					powf(cabsf(in[x*ks*ls + (ls-z)]), 2.)) ;
				
				out[x*ks*ls + z] = ave_intens ;
				out[x*ks*ls + (ls-z)] = ave_intens ;
			}
			
			for (y = 1 ; y <= kc ; ++y) {
				ave_intens = 0.5 * (
					powf(cabsf(in[x*ks*ls + y*ls]), 2.) + 
					powf(cabsf(in[x*ks*ls + (ks-y)*ls]), 2.)) ;
				
				out[x*ks*ls + y*ls] = ave_intens ;
				out[x*ks*ls + (ks-y)*ls] = ave_intens ;
			}
			
			out[x*ks*ls] = powf(cabsf(in[x*ks*ls]), 2.) ;
		}
	}
}

// Bragg projection 
// (hklvol,vol,hsize,ksize,lsize,hoffset,koffset,loffset)
// (rhkl,fhkl,exp_hkl,hkl_mag)
void proj_bragg(float *in, float *out) {
	long h, k, l, i ;
	float norm_factor = 1.f / (float)hklvol ;
	
	memset(exp_hkl, 0, hklvol*sizeof(float)) ;
	memset(out, 0, vol*sizeof(float)) ;
	
	// Extract central region and rescale axes
	for (h = 0 ; h < hsize ; ++h)
	for (k = 0 ; k < ksize ; ++k)
	for (l = 0 ; l < lsize ; ++l)
		rhkl[h*ksize*lsize + k*lsize + l] 
		 = in[(h+hoffset)*size*size + (k+koffset)*size + (l+loffset)] ;
	
	// Fourier transform to get structure factors
	fftwf_execute(forward_hkl) ;
	
	// Symmetrize hkl structure factors
//	symmetrize_intens(fhkl, exp_hkl, 0) ;
	
	// Replace with measured modulus
	for (i = 0 ; i < hklvol ; ++i) {
		if (hkl_mag[i] > 0.)
//			fhkl[i] *= hkl_mag[i] / sqrt(exp_hkl[i]) ;
			fhkl[i] *= hkl_mag[i] / cabsf(fhkl[i]) ;
		else if (hkl_mag[i] == 0.)
			fhkl[i] = 0. ;
	}
	
	// Inverse Fourier transform
	fftwf_execute(inverse_hkl) ;
	
	// Zero pad and rescale axes
	for (h = 0 ; h < hsize ; ++h)
	for (k = 0 ; k < ksize ; ++k)
	for (l = 0 ; l < lsize ; ++l)
		out[(h+hoffset)*size*size + (k+koffset)*size + (l+loffset)]
		 = crealf(rhkl[h*ksize*lsize + k*lsize + l]) * norm_factor ;
}

// Projection over continuous data 
// (vol,rdensity,fdensity,exp_intens,obs_mag)
void proj_cont(float *in, float *out) {
	long i ;
	float norm_factor = 1.f / (float)vol ;
	memset(exp_intens, 0, vol*sizeof(float)) ;
	
	// Fourier transform to get structure factors
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = in[i] ;
	fftwf_execute(forward_cont) ;
	
	// Symmetrize to get intensities to compare
//	symmetrize_intens(fdensity, exp_intens, 1) ;
	
	// Scale using measured modulus at high resolution
	for (i = 0 ; i < vol ; ++i) {
		if (obs_mag[i] > 0.)
//			fdensity[i] *= obs_mag[i] / sqrt(exp_intens[i]) ;
			fdensity[i] *= obs_mag[i] / cabsf(fdensity[i]) ;
		else if (obs_mag[i] == 0.)
			fdensity[i] = 0. ;
	}
	
	// Inverse Fourier transform
	fftwf_execute(inverse_cont) ;
	
	for (i = 0 ; i < vol ; ++i)
		out[i] = crealf(rdensity[i]) * norm_factor ;
}

// Support projection
// (vol,num_supp,support)
void proj_supp(float *in, float *out) {
	long i, pixel ;
	
	memset(out, 0, vol*sizeof(float)) ;
	
	for (i = 0 ; i < num_supp ; ++i) {
		pixel = support[i] ;
//		if (in[pixel] > 0.)
			out[pixel] = in[pixel] ;
	}
}

// Divide constraint
void proj1(float **in, float **out) {
	proj_bragg(in[0], out[0]) ;
	proj_cont(in[1], out[1]) ;
	proj_supp(in[2], out[2]) ;
	
/*	int x ;
	float t0 = 0.f, t1 = 0.f, t2 = 0.f ;
	float it0 = 0.f, it1 = 0.f, it2 = 0.f ;
	
	for (x = 0 ; x < vol ; ++x) {
		t0 += out[0][x] ;
		t1 += out[1][x] ;
		t2 += out[2][x] ;
		it0 += in[0][x] ;
		it1 += in[1][x] ;
		it2 += in[2][x] ;
	}
	
	fprintf(stderr, "itotals = %.3e %.3e %.3e\n", it0, it1, it2) ;
	fprintf(stderr, "totals = %.3e %.3e %.3e\n", t0, t1, t2) ;
*/
}

// Concur constraint 
// (vol)
void proj2(float **in, float **out) {
	long i ;
	float ave ;
	
	for (i = 0 ; i < vol ; ++i) {
		ave = (in[0][i] + in[1][i] + in[2][i]) / 3. ;
		
		out[0][i] = ave ;
		out[1][i] = ave ;
		out[2][i] = ave ;
	}
}

// Difference map with beta = 1 
// (p1,p2,r1,vol)
double diffmap(float **x) {
	long i, t ;
	float diff, change = 0. ;
	
	proj1(x, p1) ;
	
	for (t = 0 ; t < 3 ; ++t)
	for (i = 0 ; i < vol ; ++i)
		r1[t][i] = 2. * p1[t][i] - x[t][i] ;
	
	proj2(r1, p2) ;
	
	for (t = 0 ; t < 3 ; ++t)
	for (i = 0 ; i < vol ; ++i) {
		diff = p2[t][i] - p1[t][i] ;
		x[t][i] += diff ;
		change += diff*diff ;
	}
	
	return sqrt(change / (3 * vol)) ;
}
