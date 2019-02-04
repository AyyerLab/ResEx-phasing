#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "../../src/utils.h"
#include "../../src/volume.h"
#include "../../src/quat.h"

int main(int argc, char *argv[]) {
	float *input_array, *output_array ;
	double sigma ;
	long vol, size ;
	int num_div ;
	FILE *fp ; 
	char fname_out[500] ;
	struct rotation *quat ;
	struct volume_data *volume ;
	
	/*Parse arguments*/
	if (argc < 4) {
		fprintf(stderr, "Format: %s <intens_fname> <sigma_deg> <num_div>\n", argv[0]) ;
		return 1 ;
	}
	sigma = atof(argv[2])*0.017453 ; /*degrees to radians*/
	size = get_size(argv[1], sizeof(float)) ;
	vol = size*size*size ;
	num_div = atof(argv[3]) ;
	
	/* Read input*/
	input_array = malloc(vol * sizeof(float)) ;
   	fp = fopen(argv[1], "rb") ;
	fread(input_array, sizeof(float), vol, fp) ;
	fclose(fp) ;
	
	/*Sample quad angles */
	quat = malloc(sizeof(struct rotation)) ;
	quat_gen(quat, num_div, sigma) ; 
	
	/* Blur Volume */
	volume = malloc(sizeof(struct volume_data)) ;
	volume->size = size ;
	output_array = malloc(vol * sizeof(float)) ;
	volume_rotational_blur(volume, input_array, output_array, quat) ;
	
	/* Save output */
	sprintf(fname_out, "%s-%.1fdegrees.raw", remove_ext(argv[1]), atof(argv[2])) ;
	fp = fopen(fname_out, "wb") ;
 	fwrite(output_array, sizeof(float), vol, fp) ;   
	fclose(fp) ;
	
	return 0 ;
}



