/*	
	step_A.c
*/

#include "step_A.h"

void init_mpz_array(mpz_t *array, size_t s_array){
      size_t i; 

      for (i = 0; i < s_array; i++){
            mpz_init(array[i]);
      }
}

void free_mpz_array(mpz_t *array, size_t s_array){
      size_t i; 
      for (i = 0; i < s_array; i++){
            mpz_clear(array[i]);
      }
      free(array); 
      array = NULL;
}
