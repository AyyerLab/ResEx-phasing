#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	long x[3], i[3], t, s = 501, c = 250, vol = s*s*s ;
	long a[3] = {3, 4, 6} ;
	float *model ;
	FILE *fp ;
	char fname[500] ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <model_fname>\n", argv[0]) ;
		return 1 ;
	}
	
	fp = fopen(argv[1], "rb") ;
	model = malloc(vol * sizeof(float)) ;
	fread(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	// For each reciprocal lattice point
	for (i[0] = -c/a[0] ; i[0] < c/a[0] ; ++i[0])
	for (i[1] = -c/a[1] ; i[1] < c/a[1] ; ++i[1])
	for (i[2] = -c/a[2] ; i[2] < c/a[2] ; ++i[2]) {
		for (t = 0 ; t < 3 ; ++t)
			x[t] = i[t] * a[t] + c ;
		
		if (x[0] < 1 || x[0] > s-2
		 || x[1] < 1 || x[1] > s-2
		 || x[2] < 1 || x[2] > s-2)
			continue ;
		
		model[x[0]*s*s + x[1]*s + x[2]] = (
			model[(x[0]+1)*s*s + x[1]*s + x[2]] +
			model[(x[0]-1)*s*s + x[1]*s + x[2]] +
			model[x[0]*s*s + (x[1]+1)*s + x[2]] +
			model[x[0]*s*s + (x[1]-1)*s + x[2]] +
			model[x[0]*s*s + x[1]*s + (x[2]+1)] +
			model[x[0]*s*s + x[1]*s + (x[2]-1)]
			) / 6. ;
	}
	
	sprintf(fname, "%s_nb.raw", remove_ext(argv[1])) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	
	return 0 ;
}
