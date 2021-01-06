/*
	lp_var.h
*/

#ifndef LP_VAR_H
#define LP_VAR_H

#include <gmp.h>
#include "step_A.h"
#include "step_B.h"


typedef struct AQp_lp{
      /* The goal of this structure is to create a linked list whose
       * nodes store the data of an A-Q pairs (Anm1, Qn), where Qn 
       * is almost completely factorisable with the primes of the factor
       * base. The remaining cofactor is called lp (for large prime). 
       * This linked list has two properties. It is stored according
       * to the member lp. Only one AQp_lp exists for a given value of
       * lp. That AQp_lp pair is used as the pivot of the Gaussian 
       * elimination to eliminate the large prime lp if another pair
       * with the same large prime cofactor is encounterded. 
       */
      mpz_t Qn; 
      mpz_t Anm1; 
      mpz_t lp;               // The large prime. 
      mpz_t exp_vect;         // The exponent vector associated to Qn
                              // (the large prime lp is missing).
      struct AQp_lp *next;    // The next node of the linked list. 
}AQp_lp; 


int is_qn_fact_lp_var(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows, const mpz_t Qn,
                      mpz_t lp, const mpz_t *factor_base, const size_t s_fb,
                      const mpz_t pm_squared);

AQp_lp *create_AQp_lp(const mpz_t Qn, const mpz_t Anm1, const mpz_t lp, Data_exp_vect *D_exp_vect, size_t n); 

void delete_AQp_lp_list(AQp_lp **list);


AQp_lp *insert_or_eliminate_lp(AQp_lp *list, const mpz_t Qn, const mpz_t Anm1,
                              const mpz_t lp, Data_exp_vect *D_exp_vect, size_t n,
                              mpz_t *Qns, mpz_t *Ans, mpz_t *exp_vects, size_t *nb_AQp,
                              const mpz_t N, mpz_t A, mpz_t Q, mpz_t gcd, mpz_t exp_vect); 

int create_AQ_pairs_lp_var(const Params P, mpz_t *Ans, mpz_t *Qns, size_t *nb_AQp,
                            mpz_t *exp_vects, const mpz_t *factor_base, 
                            const size_t s_fb, AQp_lp **list, mpz_t fact_found); 


#endif 
