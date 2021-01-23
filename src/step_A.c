/* step_A.c */

#include "step_A.h"

void init_mpz_array(mpz_t *array, size_t s_array) {
    for (size_t i = 0; i < s_array; i++) {
        mpz_init(array[i]);
    }
}

void free_mpz_array(mpz_t *array, size_t s_array) {
    for (size_t i = 0; i < s_array; i++) {
        mpz_clear(array[i]);
    }
    free(array); 
    array = NULL;
}

void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N,
					  unsigned k) {
    /*
    Compute the factor base of size s_fb, which contains the prime 2 and
    the smallest primes p such that the Legendre symbol (kN/p) = 1 or 0. 

    param factor_base: An array of size 's_fb' already allocated but
                       not initialized.
    param s_fb: The size of the factor base wanted.
    param N: The integer to be factored.
    param k: The multiplier k used in the expansion of sqrt(kN).
    */

    mpz_t  kN;
    mpz_t  prime; 
    size_t i;

    mpz_init(kN);
    mpz_init_set_ui(prime, 2); 
    mpz_init_set_ui(factor_base[0], 2);

    mpz_mul_ui(kN, N, k);
    i = 1; 
    while (i < s_fb) {
        mpz_nextprime(prime, prime); 
        if (-1 != mpz_legendre(kN, prime)) {
            mpz_init_set(factor_base[i], prime); 
            i ++; 
        }
    }

    mpz_clears(kN, prime, NULL);
}

void init_results(Results *R) {
    /*
    Initialize a structure Results.

    param R: Pointer to the structure to be initialized.
    */

    mpz_init(R-> fact_found); 
    R-> found  = 0; 
    R-> nb_AQp = 0; 
    R-> n_last = 0; 
} 

void clear_Params_Results(Params *P, Results *R) {
    mpz_clears(P-> N, R-> fact_found, NULL); 
}
