/* step_A.h */

#ifndef STEP_A
#define STEP_A

#include <gmp.h>
#include <stdlib.h>

struct Params {
	// Integer to factor
	mpz_t	 N;
	// Multiplier k used in sqrt(kN)
	unsigned k;
	// Upper limit on n for the cont. frac. expansion of sqrt(kN)
	size_t   n_lim;
	// Number of wanted A-Q pairs where Qn is factored in factor_base
	size_t	 nb_want_AQp;}
;

typedef struct Params Params;

#define malloc_mpz_array(s_array) (mpz_t *)malloc((s_array)*sizeof(mpz_t))

void init_mpz_array(mpz_t *array, size_t s_array);
void free_mpz_array(mpz_t *array, size_t s_array);
void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N, unsigned k);

#endif
