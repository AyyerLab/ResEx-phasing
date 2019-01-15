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
