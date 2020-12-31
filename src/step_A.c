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

void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N){
      mpz_t prime; 
      size_t i;

      mpz_init_set_ui(prime, 2); 
      mpz_init_set_ui(factor_base[0], 2);

      i = 1; 
      while (i < s_fb){
            mpz_nextprime(prime, prime); 
            if ( 1 == mpz_legendre(N, prime) ){
                  mpz_init_set(factor_base[i], prime); 
                  i ++; 
            }
      }

      mpz_clear(prime); 
}
