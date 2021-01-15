/* step_A.h */

#ifndef STEP_A
#define STEP_A

#include <stdlib.h> 
#include <gmp.h>

typedef struct Params { 
     
    int      lp_var;      /* Boolean to indicate if the large prime 
                             variation is used */
	mpz_t	 N;			  /* Integer to factor */
	unsigned k;			  /* Multiplier k used in sqrt(kN) */
   	size_t   n_lim;		  /* Upper limit of the index n in the cont.
							 frac. expansion of sqrt(kN) */
    size_t   s_fb ;       /* The size of the factor base */
	size_t	 nb_want_AQp; /* Number of wanted A-Q pairs where Qn is
							 completely factored over the facor base. */
   
} Params;


typedef struct Results {

    int    found;         /* A boolean to indicate if a factor was found. */ 
    mpz_t  fact_found;    /* The factor found (if it exists). */
    size_t nb_AQp;        /* The number of A-Q pairs found with Qn 
                             factored over the factor base */
    size_t n_last;        /* The bigger index n in the expansion reached. */ 

} Results;

#define MALLOC_MPZ_ARRAY(s_array) (mpz_t *)malloc((s_array)*sizeof(mpz_t))

void init_mpz_array(mpz_t *array, size_t s_array);
void free_mpz_array(mpz_t *array, size_t s_array);
void init_factor_base(mpz_t *factor_base, size_t s_fb, const mpz_t N, unsigned k);
void init_results(Results *R);
void clear_Params_Results(Params *P, Results *R); 

#endif
