#include "utils.h"

char* remove_ext(char *fullName) {
	char *out = malloc(500 * sizeof(char)) ;
	strcpy(out,fullName) ;
	if (strrchr(out,'.') != NULL)
		*strrchr(out,'.') = 0 ;
	return out ;
}

char* extract_fname(char* fullName) {
	return 
		strrchr(fullName,'/') != NULL
			? strrchr(fullName,'/') + 1
			: fullName ;
}

long get_size(char *vol_fname, size_t type_size) {
	struct stat st ;
	
	stat(vol_fname, &st) ;
	long size = round(pow(((double) st.st_size) / (double) type_size, 1/3.)) ;
	fprintf(stderr, "Calculated volume size = %ld\n", size) ;
	return size ;
}

