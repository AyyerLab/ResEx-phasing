#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <libgen.h>

struct rotation {
	int num_rot, num_rot_p ;
	double *quat ;
	int icosahedral_flag ;
} ;

int quat_gen(struct rotation*, int) ;
int quat_parse(struct rotation*, char*) ;
void quat_divide(struct rotation*, int, int) ;
void quat_free(struct rotation*) ;

