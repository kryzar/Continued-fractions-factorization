/* step_2.h */

#ifndef STEP_2_H
#define STEP_2_H

#include <gmp.h>
#include "init_algo.h"

void gauss_elimination(mpz_t *exp_vects, mpz_t *hist_vects, size_t *lin_rel_indexes,
                       size_t *nb_lin_rel, size_t nb_AQp);

void calculate_A_Q(mpz_t A, const mpz_t *Ans, mpz_t Q, const mpz_t *Qns,
                  mpz_t hist_vect, const mpz_t N);

void find_factor(Results *R, const mpz_t *Ans, const mpz_t *Qns, mpz_t *exp_vects, 
                mpz_t *hist_vects, const mpz_t N);


#endif 
