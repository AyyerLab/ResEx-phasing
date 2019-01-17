#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <float.h>
#include <string.h>
#include <math.h>

struct ccp4_header {
	int32_t nx, ny, nz ;
	int32_t mode ;
	int32_t nxstart, nystart, nzstart ;
	int32_t mx, my, mz ;
	float xlen, ylen, zlen ;
	float alpha, beta, gamma ;
	int32_t mapc, mapr, maps ;
	float dmin, dmax, dmean ;
	int32_t ispg ;
	int32_t nsymbt ;
	char extra[100] ;
	float xorig, yorig, zorig ;
	char cmap[4] ;
	char machst[4] ;
	float rms ;
	int32_t nlabl ;
	char labels[10][80] ;
} ;

struct ccp4_map {
	struct ccp4_header header ;
	char *ext_header ;
	
	int8_t *i8_data ;
	int16_t *i16_data ;
	float *data ;
	float complex *c8_data ;
	uint16_t *u16_data ;
} ;

int parse_map(char*, struct ccp4_map*) ;
int write_map(char*, struct ccp4_map*) ;
int save_vol_as_map(char*, float*, int[3], float[3], char*, int) ;
int save_mask_as_map(char*, int8_t*, int[3], float[3], char*, int) ;
void free_map(struct ccp4_map*) ;
