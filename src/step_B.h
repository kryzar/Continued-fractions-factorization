/* step_B.h */

#ifndef STEP_B_H
#define STEP_B_H

#include <gmp.h>
#include "step_A.h"
#include "lp_var.h"

 
int is_Qn_factorisable(const Params *P, size_t *Qn_odd_pows, 
                       size_t *nb_Qn_odd_pows, const mpz_t Qn, mpz_t Qn_divided,
                       const mpz_t *factor_base, const mpz_t pm_squared); 

void init_hist_vects(mpz_t *hist_vects, const size_t nb_AQp);

void init_exp_vect(int init, mpz_t exp_vect, const size_t *Qn_odd_pows, 
                   size_t nb_Qn_odd_pows, size_t n); 

void create_AQ_pairs(const Params *P, Results *R, mpz_t *Ans, mpz_t *Qns, 
                     mpz_t *exp_vects, const mpz_t *factor_base, 
                     AQp_lp **list);

#endif
