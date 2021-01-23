/* fact.h */

#ifndef FACT_H
#define FACT_H

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <time.h> 
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"

#define S_PRIMES 10  // Size of the PRIMES array used to choose_k 
#define K_MAX    97  // Default value for the maximum value of k

void print_results(const Params *P, const Results *R);

void contfract_factor(const Params *P, Results *R);

size_t choose_s_fb(const mpz_t N); 

unsigned choose_k(const mpz_t N, unsigned k_max); 

#endif 
