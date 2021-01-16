/* fact.c */

#include "fact.h"

void print_results(const Params *P, const Results *R){

    // Print the parameters
    if (P -> lp_var) {
        printf("With the large prime variation \n"); 
    }else{
        printf("Without the large prime variation \n"); 
    }
    gmp_printf("N: %Zd \n",		  P -> N); 
    printf("k: %u \n",			  P -> k); 
    printf("n_lim: %lu \n",		  P -> n_lim); 
    printf("s_fb : %lu \n",		  P -> s_fb); 
    printf("nb_want_AQp: %lu \n", P -> nb_want_AQp);

    // Print the results
    if (R -> found) {
       gmp_printf("Factor found: %Zd \n", R -> fact_found); 
    }else{
        printf("No factor found \n"); 
    }
    printf("nb_AQp: %lu \n", R -> nb_AQp);
    printf("last_n: %lu \n", R -> n_last); 
} 

void contfract_factor(const Params *P, Results *R){

    mpz_t *factor_base; 
    mpz_t *Ans; 
    mpz_t *Qns; 
    mpz_t *exp_vects; 
    mpz_t *hist_vects; 	
    AQp_lp* list_AQp_lp = NULL; // For the large prime variation

    /**************
    * Allocations *
    **************/

    factor_base = MALLOC_MPZ_ARRAY(P -> s_fb);
    Ans = MALLOC_MPZ_ARRAY(P -> nb_want_AQp); 
    Qns = MALLOC_MPZ_ARRAY(P -> nb_want_AQp);
    exp_vects = MALLOC_MPZ_ARRAY(P -> nb_want_AQp); 
    hist_vects = MALLOC_MPZ_ARRAY(P -> nb_want_AQp);

    /************************
    *  Looking for a factor *
    *************************/

    // Initialize the factor base
    init_factor_base(factor_base, P -> s_fb, P -> N, P -> k);

    // Create the A-Q pairs
    if (P -> lp_var) {
        create_AQ_pairs_lp_var(P, R, Ans, Qns, exp_vects, factor_base, &list_AQp_lp); 
    }else{
        create_AQ_pairs(P, R, Ans, Qns, exp_vects, factor_base); 
    }

    // Try to find a factor
    if (! R->found){
        init_hist_vects(hist_vects, R->nb_AQp);
        find_factor(R, Ans, Qns, exp_vects, hist_vects, P->N);
    }

    /*******
    * Free *
    *******/

    free_mpz_array(Ans, R -> nb_AQp); 
    free_mpz_array(Qns, R -> nb_AQp); 
    free_mpz_array(exp_vects, R -> nb_AQp); 
    free_mpz_array(hist_vects, R -> nb_AQp); 
    free_mpz_array(factor_base, P -> s_fb);
    if (P -> lp_var){
        // For the large prime variation 
        delete_AQp_lp_list(&list_AQp_lp); 
    }
}
