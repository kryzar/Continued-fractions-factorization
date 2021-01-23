/* step_A.c */

#include "step_A.h"
#include <gmp.h>

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

void set_eas_params(Params *P, unsigned eas_cut, unsigned long eas_coeff) {
    /* 
    Indicate that the early abort strategy is used: set P-> eas to 1.
    Initialize and set the parameters P-> eas_cut, P-> eas_coeff and 
    P-> eas_bound_div of the early abort strategy.

    param P: A pointer to the structure containing the parameters. 
             P-> eas_bound_div is already initialized.
    param eas_cut: The index of the cut. If aes_cut = 0, take the default
                   value EAS_CUT.
    param eas_coeff: coefficient used to set P-> eas_bound_div to 
                     [sqrt(kN)]/eas_coeff. If bound = 0, take the default
                     value EAS_COEFF.
    */
    mpz_t kN; 

    mpz_init(kN);

    mpz_mul_ui(kN, P-> N, P-> k); 
    P-> eas = 1; 
    if (eas_cut) {
        P-> eas_cut = eas_cut;
    } else {
        P-> eas_coeff = EAS_CUT;  
    }
    if (eas_coeff) {
        P-> eas_coeff = eas_coeff; 
    } else {
        P-> eas_coeff = EAS_COEFF; 
    }
    mpz_sqrt(P-> eas_bound_div, kN);
    mpz_cdiv_q_ui(P-> eas_bound_div, P-> eas_bound_div, P-> eas_coeff); 

    mpz_clear(kN); 
}

void init_Params_Results(Params *P, Results *R) {
    /*
    Initialize the mpz_t (s) of the structures P and R.
    Set the other R's value to 0. Set by default P-> eas to 0. If
    set_aes_params is called, its value will be set to 1. 
    */

    mpz_inits(P-> N, P-> eas_bound_div, R-> fact_found, NULL);
    P-> eas    = 0; 
    R-> found  = 0; 
    R-> nb_AQp = 0; 
    R-> n_last = 0; 
    R-> time   = 0; 
}

void clear_Params_Results(Params *P, Results *R) {
    mpz_clears(P-> N, P-> eas_bound_div, R-> fact_found, NULL); 
}
