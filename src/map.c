#include <stdio.h>
#include <stdlib.h>

struct ccp4_map {
	int nx, ny, nz ;
	int mode ;
	int nxstart, nystart, nzstart ;
	int mx, my, mz ;
	float xlen, ylen, zlen ;
	float alpha, beta, gamma ;
	int mapc, mapr, maps ;
	float dmin, dmax, dmean ;
	int ispg ;
	int nsymbt ;
	char extra[100] ;
	float xorg, yorg, zorg ;
	char cmap[4] ;
	char machst[4] ;
	float rms ;
	int nlabl ;
	char labels[10][80] ;
} ;


