/*
	step_B.h
*/

#ifndef STEP_B_H
#define STEP_B_H

#include <gmp.h>
#include "step_A.h"


typedef struct Data_exp_vect{

      /* This structure gathers all the data needed to compute the exponent
       * vector of a given Qn, except the subscript n. */

      size_t *reduced_fb_indexes;    /* The elements of this array indicate
                                        the index of the primes of the factor
                                        base that have an odd valuation in at 
                                        least one of the Qn's factorization
                                        (with Qn completely factored or almost
                                        completely factored if the large prime 
                                        variation is used).                            */ 
      size_t  nb_reduced_fb_indexes; /* The number of indexes added to the 
                                        reduced_fb_indexes array.            */
      size_t *Qn_odd_pows;           /* It stores the indexes of the prime 
                                        factors of Qn that have an odd valuation.*/
      size_t nb_Qn_odd_pows;         // The number of such prime factors of Qn.
     
}Data_exp_vect;
 
int is_qn_factorisable(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows, const mpz_t Qn,
                       mpz_t Q_temp, const mpz_t *factor_base, const size_t s_fb);

void init_hist_vects(mpz_t *hist_vects, const size_t nb_AQp);

void init_exp_vect(const int init, mpz_t exp_vect, Data_exp_vect *D, const size_t n); 

int create_AQ_pairs(const Params P, mpz_t *Ans, mpz_t *Qns, size_t *nb_AQp, mpz_t *exp_vects, 
                    const mpz_t *factor_base, const size_t s_fb, mpz_t fact_found);

#endif
