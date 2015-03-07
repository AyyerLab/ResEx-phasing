#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define HMAX 39
#define KMAX 65
#define LMAX 88

int main(int argc, char *argv[]) {
	int h, k, l, hsize, ksize, lsize ;
	double f, p, arg ;
	double complex *hklarray ;
	char line[999] ;
	FILE *fp ;
	
	if (argc < 2) {
		fprintf(stderr, "%s <hkl_fname>\n", argv[0]) ;
		return 1 ;
	}
	
	hsize = 2*HMAX + 1 ;
	ksize = 2*KMAX + 1 ;
	lsize = 2*LMAX + 1 ;
	
	hklarray = calloc(hsize*ksize*lsize, sizeof(double complex)) ;
	
	fp = fopen(argv[1], "r") ;
	while (fgets(line, 999, fp) != NULL) {
		sscanf(line, "%d %d %d %lf %lf\n", &h, &k, &l, &p, &f) ;
		p *= M_PI / 180. ;
		
		arg = p ;
		hklarray[(HMAX+h)*ksize*lsize + (KMAX+k)*lsize + (LMAX+l)]
			= f * cexp(I*arg) ;
		hklarray[(HMAX-h)*ksize*lsize + (KMAX-k)*lsize + (LMAX-l)]
			= f * cexp(-I*arg) ;
		
		arg = p + (h%2 + k%2)*M_PI ;
		hklarray[(HMAX+h)*ksize*lsize + (KMAX-k)*lsize + (LMAX-l)]
			= f * cexp(I*arg) ;
		hklarray[(HMAX-h)*ksize*lsize + (KMAX+k)*lsize + (LMAX+l)]
			= f * cexp(-I*arg) ;
		
		arg = p + (k%2 + l%2)*M_PI ;
		hklarray[(HMAX-h)*ksize*lsize + (KMAX+k)*lsize + (LMAX-l)]
			= f * cexp(I*arg) ;
		hklarray[(HMAX+h)*ksize*lsize + (KMAX-k)*lsize + (LMAX+l)]
			= f * cexp(-I*arg) ;
		
		arg = p + (l%2 + h%2)*M_PI ;
		hklarray[(HMAX-h)*ksize*lsize + (KMAX-k)*lsize + (LMAX+l)]
			= f * cexp(I*arg) ;
		hklarray[(HMAX+h)*ksize*lsize + (KMAX+k)*lsize + (LMAX-l)]
			= f * cexp(-I*arg) ;
	}
	fclose(fp) ;
	
	fp = fopen("hkl_complex.bin", "wb") ;
	fwrite(hklarray, sizeof(double complex), hsize*ksize*lsize, fp) ;
	fclose(fp) ;
	
	free(hklarray) ;
	
	return 0 ;
}
