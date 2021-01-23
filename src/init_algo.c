/* init_algo.c */

#include "init_algo.h"
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

const unsigned PRIMES[S_PRIMES] = {3, 5, 7, 11, 13, 17, 19, 23, 29, 31};

unsigned choose_k(const mpz_t N, unsigned k_max) {
    /*
    Choose a value for k between 1 and k_max according to Morrison and
    Brillhart's method. For each k in the range which allows 3 or 5 to
    be in the factor base, determine the number of primes p_i contained
    in the array PRIMES of size S_PRIMES such that the Legendre symbol 
    (kN/pi) = 0 or 1. Pick the k allowing the largest number of such 
    primes. If multiple k are possible, pick the smallest k having the
    largest sum (1/pi), with pi accepted in the factor base.

    return: The choosen k 
    param N: Number to be factored
    param k_max: The maximum k allowed
    */
    mpz_t    kN; 
    unsigned nb_p_accepted[k_max + 1]; // nb_p_accepted[k] contains the 
                                       // number of primes such that 
                                       // (kN/p) = 0 or 1 (the case 0 is 
                                       // not used).
    double   sum_1_over_p [k_max + 1]; // sum_1_over_p[k] contains the
                                       // sum of 1/pi with pi accepted
                                       // if k is choosen (the case 0 is
                                       // not used).
    unsigned max_nb_p_accepted; 
    unsigned max_sum; 
    unsigned k_best; 

    mpz_init(kN);

    max_nb_p_accepted = 0;
    // Fill the nb_p_accepted and sum_1_over_p arrays.
    for (unsigned k = 1; k < k_max + 1; k++) {
        mpz_mul_ui(kN, N, k);
        nb_p_accepted[k] = 0; 
        sum_1_over_p[k]  = 0; 
        if (    (-1 != mpz_kronecker_ui(kN, 3)) 
             || (-1 != mpz_kronecker_ui(kN, 5)) ) {
            for (unsigned i = 0; i < S_PRIMES; i++) {
                if (-1 != mpz_kronecker_ui(kN, PRIMES[i])) {
                    nb_p_accepted[k] ++; 
                    sum_1_over_p[k] += 1/(double)PRIMES[i]; 
                }
            }
            if (nb_p_accepted[k] > max_nb_p_accepted) {
                max_nb_p_accepted = nb_p_accepted[k]; 
            }
        } 
    }

    max_sum = 0;
    // Pick k which maximizes nb_p_accepted[k]. If multiple k are 
    // possible, pick the one which maximizes sum_1_over_p[k]
    for (unsigned k = 1; k < k_max + 1; k++) {
        if (max_nb_p_accepted == nb_p_accepted[k]) {
            if ( sum_1_over_p[k] > max_sum ) {
                k_best = k; 
                max_sum = sum_1_over_p[k]; 
            }
        }
    }
    
    mpz_clear(kN); 
    return k_best; 
}

size_t choose_s_fb(const mpz_t N) {
    /*
    Choose the size of the factor base according to the parameters 
    of Morrison and Brillhart's paper.

    return: The size choosen for the factor base.
    param N: The integer to be factored.
    */
    size_t nb_digits = mpz_sizeinbase(N, 10); 

    if (nb_digits <= 20) {
        return 60; 
    }
    if (nb_digits <= 23) {
        return 150; 
    }
    if (nb_digits <= 25) {
        return 200; 
    }
    if (nb_digits <= 28) {
        return 300; 
    }
    if (nb_digits <= 30) {
        return 400; 
    }
    if (nb_digits <= 32) {
        return 450; 
    }
    if (nb_digits <= 34) {
        return 500; 
    }
    if(nb_digits <= 36) {
        return 550; 
    }
    if(nb_digits <= 38) {
        return 600; 
    }
    if(nb_digits <= 40) {
        return 650; 
    }
    return 750;     
} 
