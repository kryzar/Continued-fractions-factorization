/* init_algo.h */

#ifndef INIT_ALGO_H
#define INIT_ALGO_H

#include <stdlib.h> 
#include <gmp.h>

#define EAS_CUT   50        /* Default value to initialize eas_cut         */
#define EAS_COEFF 1000000   /* Default value used to compute eas_bound_div */
#define K_MAX     97        /* Default value for the maximum value of k    */
#define N_LIM     15000000  /* Default value for n_lim                     */
#define DELTA     15        /* To choose nb_want_AQp := s_fb + DELTA       */ 
#define S_PRIMES  10        /* Size of the PRIMES array used to choose k   */

typedef struct Params {
    /*
    Store all the parameters to apply the continued fraction method for 
    factoring integers. 
    */ 
    int       lp_var;        /* Boolean to indicate if the large prime
                                variation is used */
    int       eas;           /* Boolean to indicate if the early abort
                                strategy (eas) is used */
    unsigned  eas_cut;       /* Index of the cut if the eas is used */
    unsigned long eas_coeff; /* Coefficient to compute eas_bound_div */ 
    mpz_t     eas_bound_div; /* Bound used to test if the unfactored 
                                portion of Qn is too large to continue 
                                the factorization (sqrt(kN)/eas_coeff) */   
    mpz_t     N;             /* Integer to factor */
    unsigned  k;             /* Multiplier k used in sqrt(kN) */
    size_t    n_lim;         /* Upper limit of the index n in the cont.
                                frac. expansion of sqrt(kN) */
    size_t    s_fb ;         /* The size of the factor base */
    size_t    nb_want_AQp;   /* Number of wanted A-Q pairs such that Qn
                                is factor_base-smooth (or almost smooth,
                                with possibly a squared large prime in
                                its decomposition if the large prime
                                variation is used) */
} Params;

typedef struct Results {
    /*
    Store the factor found if it exits and some data to analyse the 
    method performance. 
    */
    mp_bitcnt_t nb_bits;       /* Number of bits of N */
    int         found;         /* A boolean to indicate if a factor was
                                  found. */ 
    mpz_t       fact_found;    /* The factor found (if it exists). */
    size_t      nb_AQp;        /* The number of A-Q pairs found such 
                                  that Qn is factor_base-smooth (or 
                                  almost smooth, with possibly a squared
                                  large prime in its decomposition if 
                                  the large prime variation is used) */
    size_t      n_last;        /* The bigger index n reached */
    float       time;          /* Execution time of the function 
                                  contfract_factor in seconds */
} Results;

#define MALLOC_MPZ_ARRAY(s_array) (mpz_t *)malloc((s_array)*sizeof(mpz_t))

void init_mpz_array(mpz_t *array, size_t s_array);
void free_mpz_array(mpz_t *array, size_t s_array);

void init_Params_Results(Params *P, Results *R); 
void clear_Params_Results(Params *P, Results *R);

void set_eas_params(Params *P, unsigned eas_cut, unsigned long eas_coeff); 

void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N, unsigned k);

unsigned choose_k(const mpz_t N, unsigned k_max); 

size_t choose_s_fb(const mpz_t N); 

#endif
