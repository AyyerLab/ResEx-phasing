#include "volume.h"

void volume_init(struct volume_data *self, long size) {
	self->size = size ;
	fprintf(stderr, "Symmetrizing with point group: %s\n", self->point_group) ;

	self->intrad = NULL ;
	self->radavg = NULL ;
	self->radcount = NULL ;
	self->obs_radavg = NULL ;
}

/* Symmetrize intensity incoherently according to given point group
 * The array is assumed to have q=0 at (0,0,0) instead of in the center of the array
 * Function is parallelized using OpenMP and co cannot be used within parallel block
*/
void volume_symmetrize_incoherent(struct volume_data *self, float complex *in, float *out, float *bg) {
	long hs, ks, ls, kc, lc ;
	long size = self->size, vol = size*size*size ;
	
	hs = size ;
	ks = size ;
	ls = size ;
	kc = ks / 2 ;
	lc = ls / 2 ;
	
	if (strcmp(self->point_group, "222") == 0) {
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
	else if (strcmp(self->point_group, "4") == 0) {
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
	else if (strcmp(self->point_group, "1") == 0) {
		long x ;
		for (x = 0 ; x < vol ; ++x)
		if (bg != NULL)
			out[x] = sqrtf(powf(cabsf(in[x]), 2.f) + powf(bg[x], 2.f)) ;
		else
			out[x] = cabsf(in[x]) ;
	}
	else {
		fprintf(stderr, "Unrecognized point group: %s\n", self->point_group) ;
	}
}

/* Symmetrize intensity according to given point group
 * The array is assumed to have q=0 in the center of the array
 * Function is parallelized using OpenMP and co cannot be used within parallel block
*/
void volume_symmetrize_centered(struct volume_data *self, float complex *in, float *out) {
	long x, y, z ;
	long size = self->size, c = size/2, vol = size*size*size ;
	
	if (strcmp(self->point_group, "222") == 0) {
		#pragma omp parallel for default(shared) schedule(static, 1)
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			out[x*size*size + y*size + z] = 0.25 * (powf(cabsf(in[x*size*size + y*size + z]), 2.f) +
													powf(cabsf(in[x*size*size + y*size + (2*c-z)]), 2.f) +
													powf(cabsf(in[x*size*size + (2*c-y)*size + z]), 2.f) +
													powf(cabsf(in[x*size*size + (2*c-y)*size + (2*c-z)]), 2.f)) ;
		}
	}
	else if (strcmp(self->point_group, "4") == 0) {
		#pragma omp parallel for default(shared) schedule(static, 1)
		for (x = 0 ; x < size ; ++x)
		for (y = 0 ; y < size ; ++y)
		for (z = 0 ; z < size ; ++z) {
			out[x*size*size + y*size + z] = 0.25 * (powf(cabsf(in[x*size*size + y*size + z]), 2.f) +
													powf(cabsf(in[x*size*size + (2*c-z)*size + y]), 2.f) +
													powf(cabsf(in[x*size*size + (2*c-y)*size + (2*c-z)]), 2.f) +
													powf(cabsf(in[x*size*size + z*size + (2*c-y)]), 2.f)) ;
		}
	}
	else if (strcmp(self->point_group, "1") == 0) {
		#pragma omp parallel for default(shared) schedule(static, 1)
		for (x = 0 ; x < vol ; ++x)
			out[x] = powf(cabsf(in[x]), 2.f) ;
	}
	else {
		fprintf(stderr, "Unrecognized point group: %s\n", self->point_group) ;
	}
}

/* Radial average initialization
 * Calculate bin for each voxel and bin occupancy
 * Note that q=0 is at (0,0,0) and not in the center of the array
*/
void volume_init_radavg(struct volume_data *self) {
	long x, y, z, size = self->size, c = size/2, vol = size*size*size ;
	long dx, dy, dz ;
	double dist ;
	
	self->intrad = malloc(vol * sizeof(int)) ;
	self->radavg = calloc(size, sizeof(float)) ;
	self->obs_radavg = calloc(size, sizeof(float)) ;
	self->radcount = calloc(size, sizeof(float)) ;
	
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = x <= c ? x : size-x ; 
		dy = y <= c ? y : size-y ; 
		dz = z <= c ? z : size-z ; 
		dist = sqrt(dx*dx + dy*dy + dz*dz) ;
		self->intrad[x*size*size + y*size + z] = (int) dist ;
		self->radcount[self->intrad[x*size*size + y*size + z]] += 1. ;
	}
	
	for (x = 0 ; x < size ; ++x)
	if (self->radcount[x] == 0.)
		self->radcount[x] = 1. ;
}

/* Radial average calculation
 * Using previously calculated bins and bin occupancies
 * Note that q=0 is at (0,0,0) and not in the center of the array
*/
void volume_radial_average(struct volume_data *self, float *in, float *out) {
	long i, size = self->size, vol = size*size*size ;
	
	memset(self->radavg, 0, size*sizeof(float)) ;
	for (i = 0 ; i < vol ; ++i)
		self->radavg[self->intrad[i]] += in[i] ;
	
	for (i = 0 ; i < size ; ++i) {
		self->radavg[i] /= self->radcount[i] ;
		if (self->radavg[i] < 0.)
			self->radavg[i] = 0. ;
	}
	
	if (out != NULL)
	for (i = 0 ; i < vol ; ++i)
		out[i] = self->radavg[self->intrad[i]] ;
}

static void gen_rot(float rot[3][3], double *q) {
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
void volume_rotational_blur(struct volume_data *self, float *in, float *out, struct rotation *quat) {
	long vol = self->size*self->size*self->size ; 
	float *weight = calloc(vol, sizeof(float)) ;
	
	#pragma omp parallel default(shared)
	{
		long r, i, j, x, y, z, vox[3] ;
		long size = self->size, c = size / 2, vol = size*size*size ;
		float fx, fy, fz, cx, cy, cz, w, val ;
		float rot_vox[3], rot[3][3] ;
		float *priv_out = calloc(vol, sizeof(float)) ;
		float *priv_weight = calloc(vol, sizeof(float)) ;
		
		#pragma omp for schedule(dynamic,1)
		for (r = 0 ; r < quat->num_rot ; ++r) {
			if (quat->quat[r*5 + 4] < 0.)
				continue ;
			
			gen_rot(rot, &(quat->quat[r*5])) ;
			
			w = quat->quat[r*5 + 4] ;
			
			
			for (vox[0] = 0 ; vox[0] < size ; ++vox[0]) /* loop trough every pixel */
			for (vox[1] = 0 ; vox[1] < size ; ++vox[1])
			for (vox[2] = 0 ; vox[2] < size ; ++vox[2]) {
				for (i = 0 ; i < 3 ; ++i) { /* find rotated pixel position */
					rot_vox[i] = c ;
					
					for (j = 0 ; j < 3 ; ++j)
						rot_vox[i] += rot[i][j] * (vox[j] - c) ;
					
					rot_vox[i] = fmod(rot_vox[i] + size, size) ;
					
				}
				
				
				
				x = rot_vox[0] ; /* integer part of rotated coordinate, confusing syntax */
				y = rot_vox[1] ;
				z = rot_vox[2] ;
				fx = rot_vox[0] - x ; /*  fractional part of cooridnate used for linear interpolation in the next block  */
				fy = rot_vox[1] - y ;
				fz = rot_vox[2] - z ;
				cx = 1. - fx ; 
				cy = 1. - fy ;
				cz = 1. - fz ;
				
				val = in[vox[0]*size*size + vox[1]*size + vox[2]] * w ;
	
	            				
				
				priv_out[x*size*size + y*size + z]                                     += cx*cy*cz*val ;/* linear interpolation */
				priv_weight[x*size*size + y*size + z]                                  += cx*cy*cz*w ;
				priv_out[x*size*size + y*size + ((z+1)%size)]                          += cx*cy*fz*val ;
				priv_weight[x*size*size + y*size + ((z+1)%size)]                       += cx*cy*fz*w ;
				priv_out[x*size*size + ((y+1)%size)*size + z]                          += cx*fy*cz*val ;
				priv_weight[x*size*size + ((y+1)%size)*size + z]                       += cx*fy*cz*w ;
				priv_out[((x+1)%size)*size*size + y*size + z]                          += fx*cy*cz*val ;
				priv_weight[((x+1)%size)*size*size + y*size + z]                       += fx*cy*cz*w ;		
				priv_out[x*size*size + ((y+1)%size)*size + ((z+1)%size)]               += cx*fy*fz*val ;
				priv_weight[x*size*size + ((y+1)%size)*size + ((z+1)%size)]            += cx*fy*fz*w ;
				priv_out[((x+1)%size)*size*size + y*size + ((z+1)%size)]               += fx*cy*fz*val ;		
				priv_weight[((x+1)%size)*size*size + y*size + ((z+1)%size)]            += fx*cy*fz*w ;
				priv_out[((x+1)%size)*size*size + ((y+1)%size)*size + z]               += fx*fy*cz*val ;
				priv_weight[((x+1)%size)*size*size + ((y+1)%size)*size + z]            += fx*fy*cz*w ;
				priv_out[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)]    += fx*fy*fz*val ;
				priv_weight[((x+1)%size)*size*size + ((y+1)%size)*size + ((z+1)%size)] += fx*fy*fz*w ;
			}
			if (omp_get_thread_num() == 0) 
				fprintf(stderr, "\rRotating %ld/%d", r, quat->num_rot) ;
		
		}
		
		if (omp_get_thread_num() == 0)
			memset(out, 0, vol*sizeof(float)) ;
		#pragma omp barrier
		
		#pragma omp critical(priv_out)
		{
			for (x = 0 ; x < vol ; ++x) {
				out[x] += priv_out[x] ;
				weight[x] += priv_weight[x] ;
			}
		}
		
		free(priv_out) ;
		free(priv_weight) ;
	}
	
	long x ;
	for (x = 0 ; x < vol ; ++x)
	if (weight[x] > 0.)
		out[x] /= weight[x] ;
	free(weight) ;
	fprintf(stderr, "\n") ;
}

void volume_free(struct volume_data *self) {
	if (self->intrad != NULL)
		free(self->intrad) ;
	if (self->radavg != NULL)
		free(self->radavg) ;
	if (self->radcount != NULL)
		free(self->radcount) ;
	if (self->obs_radavg != NULL)
		free(self->obs_radavg) ;
}

void volume_accumulate(float *current, float *sum, long num_vox) {
	long i ;
	
	for (i = 0 ; i < num_vox ; ++i)
		sum[i] += current[i] ;
}

/* Save orthogonal central slices to file
*/
void volume_dump_slices(float *vol, char *fname, long size, int is_fourier) {
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

void volume_dump_support_slices(uint8_t *vol, char *fname, long size) {
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

/*Get mode of positive values in volume
 */
float volume_positive_mode(float *model, long size) {
	long i, valbin, maxhist = 0, hist[99] ;
	float bin[99], maxval = 0. ;
	long vol = size*size*size ;
	
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


