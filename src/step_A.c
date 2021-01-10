/* step_A.c */

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

void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N, unsigned k){

      /* This function computes the factor base of size s_fb, which 
       * contains the prime 2 and the smallest primes p such that
       * the Legendre symbol (kN/p) = 1. 
       *
       * param factor_base: An array of size 's_fb' already allocated but not initialized.
       * param s_fb: The size of the factor base wanted.
       * param N: The integer to be factored
       * param k: The muliplier k used in the expansion of sqrt(kN)
       */

      mpz_t kN;
      mpz_t prime; 
      size_t i;
      
      mpz_init(kN);
      mpz_init_set_ui(prime, 2); 
      mpz_init_set_ui(factor_base[0], 2);

      mpz_mul_ui(kN, N, k);
      i = 1; 
      while (i < s_fb){
            mpz_nextprime(prime, prime); 
            if ( 1 == mpz_legendre(kN, prime) ){
                  mpz_init_set(factor_base[i], prime); 
                  i ++; 
            }
      }

      mpz_clears(kN, prime, NULL); 
}
