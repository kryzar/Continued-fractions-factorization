/* test.c */

#include "test.h"

void rand_N(mpz_t N, mp_bitcnt_t nb_bits, gmp_randstate_t state) {
    /*
    Choose randomly a value for N such that N is a product of two 
    prime number of the same size and N is approximately nb_bits bit. 
    
    param N: Is already initialized. 
    param nb_bits: Number of bits wanted (approximately) for N.
    param state: Variable to generate a pseudo-random number.
    */
       
    mpz_t prime_1; 
    mpz_t prime_2;

    mpz_inits(prime_1, prime_2, NULL); 

    nb_bits = nb_bits /2 ; 
    while (! mpz_probab_prime_p(prime_1, 25)) {
        mpz_urandomb(prime_1, state, nb_bits);
        mpz_setbit(prime_1, nb_bits); 
    }
    while (! mpz_probab_prime_p(prime_2, 25)) {
        mpz_urandomb(prime_2, state, nb_bits);
        mpz_setbit(prime_2, nb_bits); 
    }
  
    mpz_mul(N, prime_1, prime_2); 

    mpz_clears(prime_1, prime_2, NULL); 
}

int nb_bits_vs_time(const char *file_name, mp_bitcnt_t start_bits, 
                    mp_bitcnt_t end_bits, mp_bitcnt_t incr) {
    /*
    Collect in a file data expressing the computation time as a function 
    of the size of N. 

    param file_name: Name of the file to write the data.
    param start_bits: 
    param end_bits: 
    param incr: 

    */

    FILE           *file; 
    Params          P; 
    Results         R; 
    unsigned long   seed; 
    gmp_randstate_t state; 
    
    seed = time(NULL); 
    gmp_randinit_default(state); 
    gmp_randseed_ui(state, seed);
    init_Params_Results(&P, &R);

    P.k = 1; 
    P.n_lim = 3300000;
    P.lp_var = 1; 
 
    if (( file = fopen(file_name, "w")) == NULL) {
        fprintf(stderr, "\ncannot open the file %s\n", file_name); 
        return 0; 
    } 
 
    for (mp_bitcnt_t i = start_bits; i <= end_bits; i += incr ){
        rand_N(P.N, i, state);
        set_eas_params(&P, EAS_CUT, EAS_COEFF);  
        P.s_fb = choose_s_fb(P.N); 
        P.nb_want_AQp = P.s_fb + 16; 
        R.found = 0; 
    
        contfract_factor(&P, &R);
        fprintf(file, "%lu\t %f \n", R.nb_bits, R.time); 
    }
    
    fclose(file); 
 
    gmp_randclear(state); 
    clear_Params_Results(&P, &R);

    return 1; 
}
