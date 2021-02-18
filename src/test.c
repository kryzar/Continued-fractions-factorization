/* test.c */

#include "test.h"
#include "fact.h"
#include "init_algo.h"
#include <float.h>

void fact_F7(int lp_var, int eas) {
    /*
    Test the factorization method on F7 = 2^128 + 1 and print the results.
    The parameter k is choosen with the function choose_k. If the 
    large prime variation is used, s_fb is choosen with the function
    choose_s_fb. If not, we take, according to our tests, the smallest size
    of s_fb to factor F7. 

    param lp_var: 1 if the large prime variation should be used, 0 otherwise
    param eas: 1 if the early abort strategy should be used, 0 otherwise
    */
    Params  P; 
    Results R; 

    init_Params_Results(&P, &R); 

    mpz_set_str(P.N, "340282366920938463463374607431768211457", 10);
    P.k = choose_k(P.N, K_MAX); 

    if (lp_var) {
        P.lp_var      = 1;
        P.s_fb        = choose_s_fb(P.N); 
        P.nb_want_AQp = P.s_fb + DELTA; 
    } else {
        P.lp_var      = 0;
        P.s_fb        = 950; 
        P.nb_want_AQp = P.s_fb + DELTA; 
    }

    if (eas) {
        set_eas_params(&P, EAS_CUT, EAS_COEFF);
    }

    contfract_factor(&P, &R); 
    print_results(&P, &R); 
    clear_Params_Results(&P, &R); 
}

void rand_N(mpz_t N, mp_bitcnt_t nb_bits, gmp_randstate_t state) {
    /*
    Choose randomly a value for N such that N is a product of two 
    prime numbers of the same size and N is approximately nb_bits bits. 
    
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

void fact_rand_N(mp_bitcnt_t nb_bits) {
    /*
    Choose with rand_N a value for N such that N is nb_bits bits. Call
    contfract_factor to try to find a factor of N. Print the results.
    
    param nb_bits: Number of bits wanted (approximately) for N.
    */
    Params          P; 
    Results         R; 
    gmp_randstate_t state;
    unsigned long   seed; 

    init_Params_Results(&P, &R); 
    seed = time(NULL); 
    gmp_randinit_default(state); 
    gmp_randseed_ui(state, seed);
    
    rand_N(P.N, nb_bits, state); 
    P.k    = choose_k(P.N, K_MAX);
    P.s_fb = choose_s_fb(P.N);  
    P.nb_want_AQp = P.s_fb + DELTA; 

    set_eas_params(&P, EAS_CUT, EAS_COEFF); 
    contfract_factor(&P, &R); 
    print_results(&P, &R);  

    gmp_randclear(state); 
    clear_Params_Results(&P, &R);
}

void fact_N(mpz_t N) {
    /*
    Call contfract_factor to try to find a factor of N. Print the 
    results.

    param N : Number to be factored. 
    */
    Params   P; 
    Results  R; 

    init_Params_Results(&P, &R);
    
    mpz_set(P.N, N); 
    P.k    = choose_k(P.N, K_MAX);
    P.s_fb = choose_s_fb(P.N);
    P.nb_want_AQp = P.s_fb + DELTA; 

    set_eas_params(&P, EAS_CUT, EAS_COEFF); 
    contfract_factor(&P, &R); 
    print_results(&P, &R); 

    clear_Params_Results(&P, &R);
}

size_t S_FB[9][12] = { {70 , 10 , 30 , 40 , 50 , 60 , 70 , 80 , 90 , 100, 110, 120},
                       {80 , 8  , 70 , 80 , 90 , 100, 110, 120, 130, 140, 0  , 0  },
                       {90 , 8  , 100, 110, 120, 130, 140, 150, 160, 170, 0  , 0  },
                       {100, 5  , 150, 160, 170, 180, 190, 0  , 0  , 0  , 0  , 0  },
                       {110, 5  , 200, 220, 240, 260, 280, 0  , 0  , 0  , 0  , 0  },
                       {120, 5  , 280, 310, 340, 370, 400, 0  , 0  , 0  , 0  , 0  },
                       {130, 8  , 280, 310 ,340, 370, 400, 430, 460, 490, 0  , 0  },
                       {140, 6  , 400, 450 ,500, 550, 600, 650, 0  , 0  , 0  , 0  },
                       {150, 4  , 450, 500, 600, 700, 0  , 0  , 0  , 0  , 0  , 0  }}; 

/*
Array used to test, with the function test_s_fb, several values of s_fb to
factor a integer.
S_FB[*][0] is the number of bits of the integer N to factor. 
S_FB[*][2 ...] are the values of s_fb for the test.
S_FB[*][1] is the number of s_fb value to test.
*/

int test_s_fb(const char *file_name, unsigned row_S_FB, unsigned nb_tests){
    /*
    Write at the end of file_name the execution time T of the programm 
    as a function of the size of the factor base and the integer N.
    The sizes s_fb to choose are in the array S_FB. Integers N must be 
    S_FB[row_S_FB][0] bits. One line is written per integer N. A line 
    looks like : 
    
    nb bits|N|T with s_fb = S_FB[*][2] |T with s_fb = S_FB[*][3] |...|best s_fb|best time 

    param file_name: Name of the file. 
    param row_S_FB: The row of S_FB to use. 
    param nb_tests: Number of integer N to be tested. 
                    1 test <-> 1 line in the file.
    */
    FILE           *file; 
    Params          P; 
    Results         R; 
    unsigned long   seed; 
    gmp_randstate_t state;
    size_t          best_s_fb; 
    float           best_time; 

    seed = time(NULL); 
    gmp_randinit_default(state); 
    gmp_randseed_ui(state, seed);
    init_Params_Results(&P, &R);

    if (( file = fopen(file_name, "a")) == NULL) {
        fprintf(stderr, "\nCannot open the file %s\n", file_name); 
        return 0; 
    }

    for (unsigned i = 0; i < nb_tests; i++) {
        best_time = FLT_MAX; 
        rand_N(P.N, S_FB[row_S_FB][0], state);
        set_eas_params(&P, EAS_CUT, EAS_COEFF);
        P.k = choose_k(P.N, K_MAX); 

        gmp_fprintf(file, "%lu\t %Zd\t",  mpz_sizeinbase(P.N, 2) - 1, P.N);      
        
        for (unsigned j = 2; j < 2 + S_FB[row_S_FB][1]; j++) {
            P.s_fb        = S_FB[row_S_FB][j]; 
            P.nb_want_AQp = P.s_fb + DELTA; 
            R.found       = 0;  
            contfract_factor(&P, &R); 
            if (R.found) {
                fprintf(file, "%f\t", R.time);
                if (R.time < best_time) {
                    best_time = R.time;
                    best_s_fb = S_FB[row_S_FB][j]; 
                }
            } else {
                 fprintf(file, "%f\t", FLT_MAX); 
            }
        } 
        fprintf(file, "%lu\t %f\t \n", best_s_fb, best_time); 
    }
 
    fclose(file); 
    gmp_randclear(state); 
    clear_Params_Results(&P, &R);
    return 1; 
}

int test_nb_bits_vs_time(const char *file_name, mp_bitcnt_t nb_bits, 
                         unsigned nb_tests) {
    /*
    Write at the end of file_name the execution time of the programm as
    a function of the number of bits of the integer N to factor. The 
    size of the factor base is choosen with the function choose_s_fb.
    We compute the time for k = 1 (time_k_1) and for k given by the 
    function choose_k(time_k_default). One line is written per integer
    N. A line looks like: 
    
    nb bits of N| k by default |time_k_1 | time_k_default| time_k_1/time_k_default

    param file_name: Name of the file. 
    param nb_bits: Number of bits of the integer N to factor. 
    param nb_tests: Number of integer N to be tested. 
                    1 test <-> 1 line in the file.
    */
    FILE           *file; 
    Params          P; 
    Results         R; 
    unsigned long   seed; 
    gmp_randstate_t state;
    float           time_k_1; 
    float           time_k_default; 

    seed = time(NULL); 
    gmp_randinit_default(state); 
    gmp_randseed_ui(state, seed);
    init_Params_Results(&P, &R);

    if (( file = fopen(file_name, "a")) == NULL) {
        fprintf(stderr, "\nCannot open the file %s\n", file_name); 
        return 0; 
    }
    for (unsigned i = 0; i < nb_tests; i++) {
        rand_N(P.N, nb_bits, state);
        set_eas_params(&P, EAS_CUT, EAS_COEFF);
        P.s_fb        = choose_s_fb(P.N);
        P.nb_want_AQp = P.s_fb + DELTA; 
        // First test with k = 1
        R.found = 0;
        P.k     = 1; 
        contfract_factor(&P, &R);
        if (R.found) {
            time_k_1 = R.time; 
        } else {
            time_k_1 = FLT_MAX; 
        }
        // Second test with k choosen by choose_k 
        R.found = 0;
        P.k     = choose_k(P.N, K_MAX); 
        contfract_factor(&P, &R); 
        if (R.found) {
            time_k_default = R.time; 
        } else {
            time_k_default = FLT_MAX; 
        }

        fprintf(file, "%lu\t %d\t %f\t %f\t %f \n", R.nb_bits, P.k, time_k_1,
                time_k_default, time_k_1/time_k_default); 

    }
    fclose(file);  
    gmp_randclear(state); 
    clear_Params_Results(&P, &R);
    return 1; 
}
