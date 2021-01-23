/* lp_var.h */

#ifndef LP_VAR_H
#define LP_VAR_H

#include <gmp.h>
#include "init_algo.h"

typedef struct AQp_lp{
    /*
    The goal of this structure is to create a linked list whose
    nodes store the data of an A-Q pairs (Anm1, Qn), where Qn 
    is almost factor_base-smooth. The only prime factor of Qn
    which is not in the factor_base is called lp (for large prime). 
    This linked list has two properties. It is stored according
    to the member lp. Only one AQp_lp exists for a given value of
    lp. That AQp_lp pair is used, if another pair with the same
    large prime cofactor is encountered, as the pivot of the 
    Gaussian elimination to eliminate the large prime lp. 
     */
    mpz_t Qn;            /* The Qn almost factor_base-smooth */
    mpz_t Anm1;          /* The Anm1 of the A-Q pair */
    mpz_t lp;            /* The large prime factor of Qn */ 
    mpz_t exp_vect;      /* The exponent vector associated to Qn
                            (without taking into account lp) */
    struct AQp_lp *next; /* The next node of the linked list. */ 
}AQp_lp; 


AQp_lp *create_AQp_lp(const mpz_t Qn, const mpz_t Anm1, const mpz_t lp, 
                      const mpz_t exp_vect); 

void insert_or_elim_lp(AQp_lp **list, const mpz_t Qn, const mpz_t Anm1, 
                       const mpz_t lp, mpz_t exp_vect, mpz_t *Qns,
                       mpz_t *Ans, mpz_t *exp_vects, size_t *nb_AQp,
                       const mpz_t N, mpz_t A, mpz_t Q, mpz_t gcd, Results *R); 

void delete_AQp_lp_list(AQp_lp **list);

#endif
