/* fact.c */

#include "fact.h"

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
    gmp_printf("N: %Zd \n",       P-> N); 
    printf("k: %u \n",            P-> k); 
    printf("n_lim: %lu \n",       P-> n_lim); 
    printf("s_fb : %lu \n",       P-> s_fb); 
    printf("nb_want_AQp: %lu \n", P-> nb_want_AQp);

    // Print the results
    if (R-> found) {
       gmp_printf("\nFactor found: %Zd \n", R-> fact_found); 
    } else {
        printf("No factor found \n"); 
    }
    printf("nb_AQp: %lu \n", R-> nb_AQp);
    printf("last_n: %lu \n", R-> n_last); 
} 

void contfract_factor(const Params *P, Results *R) {
    /*
    Use the continued fraction method to search a non trivial factor
    of P-> N. Compute the factor base with the auxiliary function 
    init_factor_base. Create the A-Q pairs with the auxiliary function
    create_AQ_pairs. At last, search a factor with the auxiliary function
    find_factor. 
    If a factor is found, R-> found will be 1 and the factor will be in
    R-> factor.

    param P: Pointer to the structure used to store the parameters.
    param R: Pointer to the structure which will be used to store the
             results.
    */

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
