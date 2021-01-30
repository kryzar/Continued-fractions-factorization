/* test.h */

#ifndef TEST_H
#define TEST_H
 
#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include <float.h>
#include "init_algo.h"
#include "step_1.h"
#include "step_2.h"
#include "fact.h"

void fact_F7(int lp_var, int eas); 

void rand_N(mpz_t N, mp_bitcnt_t nb_bits, gmp_randstate_t state);

void fact_rand_N(mp_bitcnt_t nb_bits); 

int test_s_fb(const char *file_name, unsigned row_S_FB, unsigned nb_tests); 

int test_nb_bits_vs_time(const char *file_name, mp_bitcnt_t nb_bits,
                         unsigned nb_tests); 

#endif
