/*
	step_B.h
*/

#ifndef STEP_B_H
#define STEP_B_H

#include <gmp.h>
#include "step_A.h"


int is_qn_factorisable(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows, const mpz_t Qn,
                       mpz_t Q_temp, const mpz_t *factor_base, const size_t s_fb);

void add_hist_vect(mpz_t *hist_vects, const size_t nb_AQp); 

void add_exp_vect(mpz_t *exp_vects, size_t *reduced_fb_indexes,
                  size_t *nb_reduced_fb_indexes, const size_t *Qn_odd_pows,
                  const size_t nb_Qn_odd_pows, const size_t nb_AQp, const size_t n);

void create_AQ_pairs(const Params P, mpz_t *Ans, mpz_t *Qns, size_t *nb_AQp,
                     mpz_t *exp_vects, mpz_t *hist_vects, 
                     const mpz_t *factor_base, const unsigned int s_fb);

#endif
