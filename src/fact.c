/* fact.c */

#include "fact.h"
#include <gmp.h>

void print_results(const Params *P, const Results *R) {
    /*
    Print the parameters used to search a factor of P-> N and the
    results.

    param P: Pointer to the structure used to store the parameters.
    param R: Pointer to the structure used to store the results.
    */

    // Print the parameters
    if (P-> lp_var) {
        printf("With the large prime variation \n"); 
    } else {
        printf("Without the large prime variation \n"); 
    }
    if (P-> eas) {
        printf("With the early abort strategy (1 cut) \n"); 
    } else {
        printf("Without the early abort strategy \n"); 
    }
    gmp_printf("N: %Zd \n", P-> N);
    if (P-> eas) {
        printf("eas_cut: %u th prime \n", P-> eas_cut); 
        printf("eas_coeff: %lu\n",        P-> eas_coeff); 
    }
    printf("k: %u \n",            P-> k); 
    printf("n_lim: %lu \n",       P-> n_lim); 
    printf("s_fb : %lu \n",       P-> s_fb); 
    printf("nb_want_AQp: %lu \n", P-> nb_want_AQp);

    // Print the results
    printf("\nsize of N: %lu bits \n", R-> nb_bits); 
    if (R-> found) {
       gmp_printf("Factor found: %Zd \n", R-> fact_found); 
    } else {
        printf("No factor found \n"); 
    }
    printf("nb_AQp: %lu \n", R-> nb_AQp);
    printf("last_n: %lu \n", R-> n_last);
    printf("time: %f \n\n",  R-> time); 
}

void contfract_factor(const Params *P, Results *R) {
    /*
    Use the continued fraction method to search a non trivial factor
    of P-> N. Compute the factor base with the auxiliary function 
    init_factor_base. Create the A-Q pairs with the auxiliary function
    create_AQ_pairs. At last, search a factor with the auxiliary function
    find_factor. 
    If a factor is found, R-> found will be 1 and the factor will be in
    R-> factor. Set the exectution time in R-> time and the size of N in
    R-> nb_bits.

    param P: Pointer to the structure used to store the parameters.
    param R: Pointer to the structure which will be used to store the
             results.
    */
    clock_t t_start = clock();
    clock_t t_end; 
    float   time; 
    mpz_t  *factor_base; 
    mpz_t  *Ans; 
    mpz_t  *Qns; 
    mpz_t  *exp_vects; 
    mpz_t  *hist_vects; 	
    AQp_lp *list_AQp_lp = NULL; // For the large prime variation

    /**************
    * Allocations *
    **************/

    factor_base = MALLOC_MPZ_ARRAY(P-> s_fb);
    Ans = MALLOC_MPZ_ARRAY(P-> nb_want_AQp); 
    Qns = MALLOC_MPZ_ARRAY(P-> nb_want_AQp);
    exp_vects = MALLOC_MPZ_ARRAY(P-> nb_want_AQp); 
    hist_vects = MALLOC_MPZ_ARRAY(P-> nb_want_AQp);

    /************************
    *  Looking for a factor *
    *************************/

    init_factor_base(factor_base, P-> s_fb, P-> N, P-> k);
    create_AQ_pairs(P, R, Ans, Qns, exp_vects, factor_base, &list_AQp_lp); 
   
    // Try to find a factor
    if (! R-> found) {
        init_hist_vects(hist_vects, R-> nb_AQp);
        find_factor(R, Ans, Qns, exp_vects, hist_vects, P-> N);
    }

    /*******
    * Free *
    *******/

    free_mpz_array(Ans, R-> nb_AQp); 
    free_mpz_array(Qns, R-> nb_AQp); 
    free_mpz_array(exp_vects, R-> nb_AQp); 
    free_mpz_array(hist_vects, R-> nb_AQp); 
    free_mpz_array(factor_base, P -> s_fb);
    delete_AQp_lp_list(&list_AQp_lp);

    t_end = clock(); 
    R-> time = (float) (t_end - t_start)/CLOCKS_PER_SEC;
    R-> nb_bits = mpz_sizeinbase(P-> N, 2) - 1; 
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
