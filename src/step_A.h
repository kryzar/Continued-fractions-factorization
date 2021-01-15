/* step_A.h */

#ifndef STEP_A
#define STEP_A

#include <gmp.h>
#include <stdlib.h>

typedef struct Params {
	mpz_t	 N;			  /* Integer to factor. */
	unsigned k;			  /* Multiplier k used in sqrt(kN). */
	size_t   n_lim;		  /* Upper limit of the index n in the cont.
							 frac. expansion of sqrt(kN). */
	size_t	 nb_want_AQp; /* Number of wanted A-Q pairs where Qn is
							 completely factored over factor_base. */
} Params;

#define MALLOC_MPZ_ARRAY(s_array) (mpz_t *)malloc((s_array)*sizeof(mpz_t))

void init_mpz_array(mpz_t *array, size_t s_array);
void free_mpz_array(mpz_t *array, size_t s_array);
void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N, unsigned k);

#endif
