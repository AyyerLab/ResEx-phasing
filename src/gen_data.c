#include "brcont.h"

int main(int argc, char *argv[]) {
	long i, x, y, z, h, k, l ;
	long dx, dy, dz, ch, ck, cl, center1, maxc ;
	float dist, qmax, qmin ;
	float *contintens, *hklintens ;
	FILE *fp ;
	
	if (argc < 3) {
		fprintf(stderr, "Format: %s <qmin> <qmax>\n", argv[0]) ;
		return 1 ;
	}
	qmin = atof(argv[1]) ;
	qmax = atof(argv[2]) ;
	
	// Setup arrays and parse files
	if (setup_gen())
		return 2 ;
	
	contintens = malloc(vol * sizeof(float)) ;
	hklintens = malloc(hklvol * sizeof(float)) ;
	
	center1 = size / 2 + 1 ;
	ch = hsize / 2 + 1 ;
	ck = ksize / 2 + 1 ;
	cl = lsize / 2 + 1 ;
	maxc = ch > ck ? ch : ck ;
	maxc = maxc > cl ? maxc : cl ;
	
	// Calculate amplitudes
	for (i = 0 ; i < vol ; ++i)
		rdensity[i] = iterate[i] ;
	
	fftwf_execute(forward_cont) ;
	
	// Symmetrize intensity
	symmetrize_intens(fdensity, exp_intens, 1) ;
	
	// Apply qmin limit
	for (x = 0 ; x < size ; ++x)
	for (y = 0 ; y < size ; ++y)
	for (z = 0 ; z < size ; ++z) {
		dx = (x + center1) % size  - center1 ;
		dy = (y + center1) % size  - center1 ;
		dz = (z + center1) % size  - center1 ;
		
		dist = sqrt(dx*dx + dy*dy + dz*dz) / center1 ;
		
		if (dist < qmin)
			exp_intens[x*size*size + y*size + z] = -1.f ;
		if (dist > size / 2.)
			exp_intens[x*size*size + y*size + z] = 0.f ;
		
		contintens[((x+size/2)%size)*size*size + ((y+size/2)%size)*size + ((z+size/2)%size)]
		 = exp_intens[x*size*size + y*size + z] ;
	}
	
	// Write to file
	fp = fopen("data/ps2contintens_427.raw", "wb") ;
	fwrite(contintens, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// Calculate hkl object
	for (h = 0 ; h < hsize ; ++h)
	for (k = 0 ; k < ksize ; ++k)
	for (l = 0 ; l < lsize ; ++l)
		rhkl[h*ksize*lsize + k*lsize + l] 
		 = iterate[(h+hoffset)*size*size + (k+koffset)*size + (l+loffset)] ;
	
	// Calculate hkl amplitudes
	fftwf_execute(forward_hkl) ;
	
	// Symmetrize hkl intensity
	symmetrize_intens(fhkl, exp_hkl, 0) ;
	
	// Apply qmax limit
	for (h = 0 ; h < hsize ; ++h)
	for (k = 0 ; k < ksize ; ++k)
	for (l = 0 ; l < lsize ; ++l) {
		dx = (h + ch) % hsize  - ch ;
		dy = (k + ck) % ksize  - ck ;
		dz = (l + cl) % lsize  - cl ;
		
		dist = sqrt(dx*dx + dy*dy + dz*dz) / maxc ;
		
		if (dist > qmax)
			exp_hkl[h*ksize*lsize + k*lsize + l] = -1.f ;
		
		hklintens[((h+hsize/2)%hsize)*ksize*lsize + ((k+ksize/2)%ksize)*lsize + ((l+lsize/2)%lsize)]
		 = exp_hkl[h*ksize*lsize + k*lsize + l] ;
	}
	
	// Write to file
	fp = fopen("data/ps2hklintens_427.raw", "wb") ;
	fwrite(hklintens, sizeof(float), hklvol, fp) ;
	fclose(fp) ;
	
	return 0 ;
}
