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
