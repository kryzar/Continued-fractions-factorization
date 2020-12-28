/* step_C.h */

#ifndef STEP_C_H
#define STEP_C_H 

#include <gmp.h>

void gauss_elimination(mpz_t *exp_vects, mpz_t *hist_vects, size_t *lin_rel_indexes,
                       size_t *nb_lin_rel, const size_t nb_AQp);

void find_A(mpz_t A, const mpz_t *Ans, mpz_t hist_vect, const mpz_t N); 

#endif 
