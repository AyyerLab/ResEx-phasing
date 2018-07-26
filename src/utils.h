#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

char *remove_ext(char*) ;
char *extract_fname(char*) ;
long get_size(char*, size_t) ;
