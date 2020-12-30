/*	
	step_A.c
*/

#include "step_A.h"

mpz_t * mpz_alloc_array(size_t s_array){
      return (mpz_t *)malloc(s_array*sizeof(mpz_t)); 
}

void mpz_init_array(mpz_t *array, size_t s_array){
      size_t i; 

      for (i = 0; i < s_array; i++){
            mpz_init(array[i]);
      }
}

void mpz_free_array(mpz_t *array, size_t s_array){
      size_t i; 
      for (i = 0; i < s_array; i++){
            mpz_clear(array[i]);
      }
      free(array); 
      array = NULL;
}
