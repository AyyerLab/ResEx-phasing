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
	return round(pow(((double) st.st_size) / type_size, 1/3.)) ;
}

