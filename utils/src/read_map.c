#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

int main(int argc, char *argv[]) {
	int nx, ny, nz, mode, mapx, mapy, mapz ;
	long psize, pvol, shx, shy, shz, x, y, z ;
	float lx, ly, lz, ax, ay, az ;
	float px, py, pz, qscale ;
	float *model, *padmodel ;
	FILE *fp ;
	char fname[999] ;
	
	if (argc < 2) {
		fprintf(stderr, "Format: %s <map_fname>\n", argv[0]) ;
		fprintf(stderr, "Optional argument: qscale (default 2200)\n") ;
		return 1 ;
	}

	qscale = 2200. ;
	if (argc == 3)
		qscale = atof(argv[2]) ;
	
	fp = fopen(argv[1], "rb") ;
	
	// Grid size
	fread(&nx, sizeof(int), 1, fp) ;
	fread(&ny, sizeof(int), 1, fp) ;
	fread(&nz, sizeof(int), 1, fp) ;
	fprintf(stderr, "Volume = (%d, %d, %d)\n", nx, ny, nz) ;
	
	// Mode
	fread(&mode, sizeof(int), 1, fp) ;
	fprintf(stderr, "Mode: %d (2 = float)\n", mode) ;
	
	fseek(fp, 12, SEEK_CUR) ;
	
	// Grid size
	int mx, my, mz ;
	fread(&mx, sizeof(int), 1, fp) ;
	fread(&my, sizeof(int), 1, fp) ;
	fread(&mz, sizeof(int), 1, fp) ;
	fprintf(stderr, "Grid size = (%d, %d, %d)\n", mx, my, mz) ;
	
	// Cell size
	fread(&lx, sizeof(float), 1, fp) ;
	fread(&ly, sizeof(float), 1, fp) ;
	fread(&lz, sizeof(float), 1, fp) ;
	
	// Cell angles
	fread(&ax, sizeof(float), 1, fp) ;
	fread(&ay, sizeof(float), 1, fp) ;
	fread(&az, sizeof(float), 1, fp) ;
	fprintf(stderr, "Cell parameters: (%.3f, %.3f, %.3f) with angles (%.1f, %.1f, %.1f)\n",
	        lx, ly, lz, ax, ay, az) ;
	
	// Mapping of axes
	fread(&mapx, sizeof(int), 1, fp) ;
	fread(&mapy, sizeof(int), 1, fp) ;
	fread(&mapz, sizeof(int), 1, fp) ;
	fprintf(stderr, "xyz -> abc mapping: (%d, %d, %d)\n", mapx, mapy, mapz) ;
	
	fseek(fp, 948, SEEK_CUR) ;
	
	model = malloc(nx*ny*nz* sizeof(float)) ;
	fread(model, sizeof(float), nx*ny*nz, fp) ;
	fclose(fp) ;
	
	sprintf(fname, "data/%s-map.raw", remove_ext(extract_fname(argv[1]))) ;
	fprintf(stderr, "Saving model to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(model, sizeof(float), nx*ny*nz, fp) ;
	fclose(fp) ;
	
	px = qscale*mx/3./lx ;
	py = qscale*my/3./ly ;
	pz = qscale*mz/3./lz ;
	px = px > py ? px : py ;
	px = px > pz ? px : pz ;
	psize = (int) px + 3 ;
	psize = psize%2 == 0 ? psize + 1 : psize ;
	pvol = psize*psize*psize ;
	
	fprintf(stderr, "Padded volume size = %ld\n", psize) ;
	px = qscale*mx/3./lx ;
	fprintf(stderr, "Ideal pad sizes = (%.3f, %.3f, %.3f)\n", pz, py, px) ;
	fprintf(stderr, "Stretch factors = (%.5f, %.5f, %.5f)\n", psize/pz, psize/py, psize/px) ;
	
	shx = (psize - nx) / 2 - 15 ;
	shy = (psize - ny) / 2 - 15 ;
	shz = (psize - nz) / 2 - 15 ;
	
	padmodel = calloc(pvol, sizeof(float)) ;
	if (mapx == 1) {
		for (x = 0 ; x < nx ; ++x)
		for (y = 0 ; y < ny ; ++y)
		for (z = 0 ; z < nz ; ++z)
			padmodel[(z+shz)*psize*psize + (y+shy)*psize + (x+shx)]
			 = model[x*ny*nz + y*nz + z] ;
	}
	else {
		for (z = 0 ; z < nz ; ++z)
		for (y = 0 ; y < ny ; ++y)
		for (x = 0 ; x < nx ; ++x)
			padmodel[(x+shx)*psize*psize + (y+shy)*psize + (z+shz)]
			 = model[((z+nz/2)%nz)*ny*nx + ((y+ny/2)%ny)*nx + ((x+nx/2)%nx)] ;
	}
	
	sprintf(fname, "data/%s-%ld.raw", remove_ext(extract_fname(argv[1])), psize) ;
	fprintf(stderr, "Saving padded model to %s\n", fname) ;
	fp = fopen(fname, "wb") ;
	fwrite(padmodel, sizeof(float), pvol, fp) ;
	fclose(fp) ;
	
	free(model) ;
	free(padmodel) ;
	
	return 0 ;
}
