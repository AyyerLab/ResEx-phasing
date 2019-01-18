#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include "quat.h"
#include "map.h"

struct volume_data {
	long size ;
	char point_group[1024] ;
	int *intrad ;
	float *radavg, *radcount, *obs_radavg ;
} ;

void volume_init(struct volume_data*, long) ;
void volume_symmetrize_incoherent(struct volume_data*, float complex*, float*, float*) ;
void volume_symmetrize_centered(struct volume_data*, float complex*, float*) ;
void volume_init_radavg(struct volume_data*) ;
void volume_radial_average(struct volume_data*, float*, float*) ;
void volume_rotational_blur(struct volume_data*, float*, float*, struct rotation*) ;
void volume_free(struct volume_data*) ;

void volume_accumulate(float*, float*, long) ;
void volume_dump_slices(float*, char*, long, int, char*) ;
void volume_dump_support_slices(int8_t*, char*, long, char*) ;
float volume_positive_mode(float*, long) ;
