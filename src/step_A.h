/*
	step_A.h
*/

#ifndef STEP_A
#define STEP_A

#include <gmp.h>
#include <stdlib.h>

struct Params {
	mpz_t	 N;			// Integer to factor
	unsigned k;	            // Integer k used in sqrt(kN)
	size_t n_lim;		// Size of the cont. frac. expansion of sqrt(kN)
	size_t nb_want_AQp;     // Number of wanted A-Q pairs (A_{n-1}, Qn) with
                              // Qn factorisable with the primes of the factor base
};

// je sais pas pour le type de n_lim 

typedef struct Params Params;

#define malloc_mpz_array(s_array) (mpz_t *)malloc((s_array)*sizeof(mpz_t))

void init_mpz_array(mpz_t *array, size_t s_array);
void free_mpz_array(mpz_t *array, size_t s_array);
void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N);

#endif
