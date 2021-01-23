/* test.h */

#ifndef TEST_H
#define TEST_H
 
#include <time.h>
#include <stdio.h>
#include <gmp.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"
#include "fact.h"

void rand_N(mpz_t N, mp_bitcnt_t nb_bits, gmp_randstate_t state); 

int nb_bits_vs_time(const char *file_name, mp_bitcnt_t start_bits, 
                    mp_bitcnt_t end_bits, mp_bitcnt_t incr); 

#endif
