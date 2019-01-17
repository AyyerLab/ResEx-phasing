#include "map.h"

int parse_map(char *fname, struct ccp4_map *map) {
	int nx, ny,nz, mode ;
	long vol ;
	FILE *fp ;
	
	fp = fopen(fname, "rb") ;
	if (fp == NULL) {
		fprintf(stderr, "Cannot open file %s\n", fname) ;
		return 1 ;
	}
	fread(&(map->header), sizeof(struct ccp4_header), 1, fp) ;
	fseek(fp, map->header.nsymbt, SEEK_CUR) ;
	
	mode = map->header.mode ;
	nx = map->header.nx ;
	ny = map->header.ny ;
	nz = map->header.nz ;
	vol = nx*ny*nz ;
	
	if (mode == 0) {
		fprintf(stderr, "Parsing (%d, %d, %d) int8_t volume\n", nx, ny, nz) ;
		map->i8_data = malloc(vol * sizeof(int8_t)) ;
		fread(map->i8_data, sizeof(int8_t), vol, fp) ;
	}
	else if (mode == 1) {
		fprintf(stderr, "Parsing (%d, %d, %d) int16_t volume\n", nx, ny, nz) ;
		map->i16_data = malloc(vol * sizeof(int16_t)) ;
		fread(map->i16_data, sizeof(int16_t), vol, fp) ;
	}
	else if (mode == 2) {
		fprintf(stderr, "Parsing (%d, %d, %d) float volume\n", nx, ny, nz) ;
		map->data = malloc(vol * sizeof(float)) ;
		fread(map->data, sizeof(float), vol, fp) ;
	}
	else if (mode == 4) {
		fprintf(stderr, "Parsing (%d, %d, %d) float complex volume\n", nx, ny, nz) ;
		map->c8_data = malloc(vol * sizeof(float complex)) ;
		fread(map->c8_data, sizeof(float complex), vol, fp) ;
	}
	else if (mode == 6) {
		fprintf(stderr, "Parsing (%d, %d, %d) uint16_t volume\n", nx, ny, nz) ;
		map->u16_data = malloc(vol * sizeof(uint16_t)) ;
		fread(map->u16_data, sizeof(uint16_t), vol, fp) ;
	}
	fclose(fp) ;
	
	return 0 ;
}

int write_map(char *fname, struct ccp4_map *map) {
	int nx, ny, nz, mode ;
	long vol ;
	FILE *fp ;
	
	fp = fopen(fname, "wb") ;
	if (fp == NULL) {
		fprintf(stderr, "Cannot open file %s\n", fname) ;
		return 1 ;
	}
	mode = map->header.mode ;
	nx = map->header.nx ;
	ny = map->header.ny ;
	nz = map->header.nz ;
	vol = nx*ny*nz ;
	
	fwrite(&(map->header), sizeof(struct ccp4_header), 1, fp) ;
	if (map->header.nsymbt > 0)
		fwrite(&(map->ext_header), sizeof(char), map->header.nsymbt, fp) ;
	
	if (mode == 0)
		fwrite(map->i8_data, sizeof(int8_t), vol, fp) ;
	else if (mode == 1)
		fwrite(map->i16_data, sizeof(int16_t), vol, fp) ;
	else if (mode == 2)
		fwrite(map->data, sizeof(float), vol, fp) ;
	else if (mode == 4)
		fwrite(map->c8_data, sizeof(float complex), vol, fp) ;
	else if (mode == 6)
		fwrite(map->u16_data, sizeof(uint16_t), vol, fp) ;
	fclose(fp) ;
	
	return 0 ;
}

int save_vol_as_map(char *fname, float *vol, int size[3], float vox_size[3], char *label, int flipped) {
	int x, nx, ny, nz ;
	float min = FLT_MAX, max = -FLT_MAX, mean = 0, rms = 0 ;
	struct ccp4_map map = {0} ;
	
	nx = size[0] ;
	ny = size[1] ;
	nz = size[2] ;
	
	for (x = 0 ; x < nx*ny*nz ; ++x) {
		mean += vol[x] ;
		rms += vol[x] * vol[x] ;
		if (vol[x] > max)
			max = vol[x] ;
		if (vol[x] < min)
			min = vol[x] ;
	}
	mean /= nx*ny*nz ;
	rms = sqrtf(rms / (nx*ny*nz) - mean*mean) ;
	
	map.header.nx = nx ;
	map.header.ny = ny ;
	map.header.nz = nz ;
	map.header.mode = 2 ;
	map.header.nxstart = 0 ;
	map.header.nystart = 0 ;
	map.header.nzstart = 0 ;
	if (flipped) {
		map.header.mx = nz ;
		map.header.my = ny ;
		map.header.mz = nx ;
		map.header.xlen = nz * vox_size[2] ;
		map.header.ylen = ny * vox_size[1] ;
		map.header.zlen = nx * vox_size[0] ;
		map.header.mapc = 3 ;
		map.header.mapr = 2 ;
		map.header.maps = 1 ;
	}
	else {
		map.header.mx = nx ;
		map.header.my = ny ;
		map.header.mz = nz ;
		map.header.xlen = nx * vox_size[0] ;
		map.header.ylen = ny * vox_size[1] ;
		map.header.zlen = nz * vox_size[2] ;
		map.header.mapc = 1 ;
		map.header.mapr = 2 ;
		map.header.maps = 3 ;
	}
	map.header.alpha = 90.f ;
	map.header.beta = 90.f ;
	map.header.gamma = 90.f ;
	map.header.dmax = max ;
	map.header.dmin = min ;
	map.header.dmean = mean ;
	map.header.ispg = 1 ;
	map.header.nsymbt = 0 ;
	memset(map.header.extra, 0, 100) ;
	map.header.xorig = 0.f ;
	map.header.yorig = 0.f ;
	map.header.zorig = 0.f ;
	strcpy(map.header.cmap, "MAP") ;
	map.header.cmap[3] = ' ' ;
	map.header.machst[0] = 0x44 ;
	map.header.machst[1] = 0x41 ; // Other two bytes zeroed
	map.header.rms = rms ;
	map.header.nlabl = 1 ;
	x = strlen(label) > 800 ? 800 : strlen(label) ;
	memcpy(map.header.labels[0], label, x) ;
	map.data = vol ;

	if (write_map(fname, &map))
		return 1 ;
	// Not freeing map because data may be used later
	
	return 0 ;
}

int save_mask_as_map(char *fname, int8_t *mask, int size[3], float vox_size[3], char *label, int flipped) {
	int x, nx, ny, nz ;
	long sum = 0 ;
	float min = FLT_MAX, max = -FLT_MAX, mean = 0, rms = 0 ;
	struct ccp4_map map = {0} ;
	
	nx = size[0] ;
	ny = size[1] ;
	nz = size[2] ;
	
	for (x = 0 ; x < nx*ny*nz ; ++x) {
		mean += mask[x] ;
		rms += mask[x] * mask[x] ;
		if (mask[x] > max)
			max = mask[x] ;
		if (mask[x] < min)
			min = mask[x] ;
		sum += mask[x] ;
	}
	mean /= nx*ny*nz ;
	rms = sqrtf(rms / (nx*ny*nz) - mean*mean) ;
	fprintf(stderr, "%ld voxels in mask\n", sum) ;
	
	map.header.nx = nx ;
	map.header.ny = ny ;
	map.header.nz = nz ;
	map.header.mode = 0 ;
	map.header.nxstart = 0 ;
	map.header.nystart = 0 ;
	map.header.nzstart = 0 ;
	if (flipped) {
		map.header.mx = nz ;
		map.header.my = ny ;
		map.header.mz = nx ;
		map.header.xlen = nz * vox_size[2] ;
		map.header.ylen = ny * vox_size[1] ;
		map.header.zlen = nx * vox_size[0] ;
		map.header.mapc = 3 ;
		map.header.mapr = 2 ;
		map.header.maps = 1 ;
	}
	else {
		map.header.mx = nx ;
		map.header.my = ny ;
		map.header.mz = nz ;
		map.header.xlen = nx * vox_size[0] ;
		map.header.ylen = ny * vox_size[1] ;
		map.header.zlen = nz * vox_size[2] ;
		map.header.mapc = 1 ;
		map.header.mapr = 2 ;
		map.header.maps = 3 ;
	}
	map.header.alpha = 90.f ;
	map.header.beta = 90.f ;
	map.header.gamma = 90.f ;
	map.header.dmax = max ;
	map.header.dmin = min ;
	map.header.dmean = mean ;
	map.header.ispg = 1 ;
	map.header.nsymbt = 0 ;
	memset(map.header.extra, 0, 100) ;
	map.header.xorig = 0.f ;
	map.header.yorig = 0.f ;
	map.header.zorig = 0.f ;
	strcpy(map.header.cmap, "MAP") ;
	map.header.cmap[3] = ' ' ;
	map.header.machst[0] = 0x44 ;
	map.header.machst[1] = 0x41 ; // Other two bytes zeroed
	map.header.rms = rms ;
	map.header.nlabl = 1 ;
	x = strlen(label) > 800 ? 800 : strlen(label) ;
	memcpy(map.header.labels[0], label, x) ;
	map.i8_data = mask ;

	if (write_map(fname, &map))
		return 1 ;
	// Not free'ing map because data may be needed later
	
	return 0 ;
}

void free_map(struct ccp4_map *map) {
	if (map->ext_header != NULL)
		free(map->ext_header) ;
	if (map->i8_data != NULL)
		free(map->i8_data) ;
	if (map->i16_data != NULL)
		free(map->i16_data) ;
	if (map->data != NULL)
		free(map->data) ;
	if (map->c8_data != NULL)
		free(map->c8_data) ;
	if (map->u16_data != NULL)
		free(map->u16_data) ;
}
