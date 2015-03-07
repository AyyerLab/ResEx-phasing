#include "brcont.h"

void init_model(float *model) {
	long i ;
	struct timeval t1 ;
	const gsl_rng_type *T ;
	gsl_rng *r ;
	
	gsl_rng_env_setup() ;
	T = gsl_rng_default ;
	gettimeofday(&t1, NULL) ;
	r = gsl_rng_alloc(T) ;
	gsl_rng_set(r, t1.tv_sec + t1.tv_usec) ;
	
	memset(model, 0, vol*sizeof(float)) ;
	
	for (i = 0 ; i < num_supp ; ++i)
		model[support[i]] = gsl_rng_uniform(r) ;
	
	gsl_rng_free(r) ;
}

void average_model(float *current, float *sum) {
	long i ;
	
	for (i = 0 ; i < vol ; ++i)
		sum[i] += current[i] ;
}

void gen_prtf(float *model) {
	long x, y, z, bin, num_bins = 50 ;
	long dx, dy, dz, center1 = size / 2 + 1 ;
	float *contrast = calloc(num_bins, sizeof(float)) ;
	long *bin_count = calloc(num_bins, sizeof(long)) ;
	FILE *fp ;
	
	// Continuous part
	for (x = 0 ; x < vol ; ++x)
		rdensity[x] = model[x] ;
	
	fftwf_execute(forward_cont) ;
	
	symmetrize_incoherent(fdensity, exp_intens) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = (x + center1) % size  - center1 ;
		dy = (y + center1) % size  - center1 ;
		dz = (z + center1) % size  - center1 ;
		
		bin = sqrt(dx*dx + dy*dy + dz*dz) / center1 * num_bins + 0.5 ;
		
		if (bin < num_bins && obs_mag[x*size*size + y*size + z] > 0.) {
//			contrast[bin] += cabs(fdensity[x*size*size + y*size + z]) / 
			contrast[bin] += sqrt(exp_intens[x*size*size + y*size + z]) / 
			                 obs_mag[x*size*size + y*size + z] ;
			bin_count[bin]++ ;
		}
	}
	
	fp = fopen("prtf.dat", "w") ;
	for (bin = 0 ; bin < num_bins ; ++bin)
	if (bin_count[bin] > 0)
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, contrast[bin]/bin_count[bin]) ;
	else 
		fprintf(fp, "%.4f\t%.6f\n", (bin + 1.) / num_bins, 0.) ;
	fclose(fp) ;
	
	free(contrast) ;
	free(bin_count) ;
}

// Symmetrize intensity incoherently according to P(2_1)(2_1)(2_1) space group
// (size)
void symmetrize_incoherent(fftwf_complex *in, float *out) {
	long hs, ks, ls, kc, lc ;
	
	hs = size ;
	ks = size ;
	ls = size ;
	kc = ks / 2 ;
	lc = ls / 2 ;
	
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

void gen_rot(float rot[3][3], int r) {
	double q0, q1, q2, q3, q01, q02, q03, q11, q12, q13, q22, q23, q33 ;
	
	q0 = quat[r*5 + 0] ;
	q1 = quat[r*5 + 1] ;
	q2 = quat[r*5 + 2] ;
	q3 = quat[r*5 + 3] ;
	
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

void blur_intens(float *in, float *out) {
	#pragma omp parallel default(shared)
	{
		long r, i, j, x, y, z, vox[3], c = size / 2 ;
		float fx, fy, fz, cx, cy, cz, w, distsq, factor ;
		float rot_vox[3], rot[3][3] ;
		float *priv_out = calloc(vol, sizeof(float)) ;
		
		#pragma omp for schedule(static,1)
		for (r = 0 ; r < num_rot ; ++r) {
			gen_rot(rot, r) ;
			
			for (vox[0] = -c ; vox[0] < size-c ; ++vox[0])
			for (vox[1] = -c ; vox[1] < size-c ; ++vox[1])
			for (vox[2] = -c ; vox[2] < size-c ; ++vox[2]) {
				distsq = vox[0]*vox[0] + vox[1]*vox[1] + vox[2]*vox[2] ;
				factor = 1. ;
				
				if (distsq < 100.*100.) {
					priv_out[((vox[0]+size)%size)*size*size + ((vox[1]+size)%size)*size + ((vox[2]+size)%size)]
					  += in[((vox[0]+size)%size)*size*size + ((vox[1]+size)%size)*size + ((vox[2]+size)%size)] *
					  	quat[r*5 + 4] ;
					continue ;
				}
				else if (distsq < 133.*133.) {
					factor = (sqrt(distsq) - 100.f) / 33.f ;
					priv_out[((vox[0]+size)%size)*size*size + ((vox[1]+size)%size)*size + ((vox[2]+size)%size)]
					  += in[((vox[0]+size)%size)*size*size + ((vox[1]+size)%size)*size + ((vox[2]+size)%size)] *
					    (1.f - factor) *
					    quat[r*5 + 4] ;
				}
				
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
				    quat[r*5 + 4] *
				    factor ;
				
				priv_out[x*size*size + y*size + z] += cx*cy*cz*w ;
				priv_out[x*size*size + y*size + ((z+1)%size)] += cx*cy*fz*w ;
				priv_out[x*size*size + ((y+1)%size)*size + z] += cx*fy*cz*w ;
				priv_out[x*size*size + ((y+1)%size)*size + ((z+1)%size)] += cx*fy*fz*w ;
				priv_out[((x+1)%size)*size*size + y*size + z] += fx*cy*cz*w ;
				priv_out[((x+1)%size)*size*size + y*size + ((z+1)%size)] += fx*cy*fz*w ;
				priv_out[((x+1)%size)*size*size + ((y+1)%size)*size + z] += fx*fy*cz*w ;
				priv_out[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += fx*fy*fz*w ;
			}
		}
		
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
