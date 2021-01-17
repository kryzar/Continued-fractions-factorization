/* fact.h */

#ifndef FACT_H
#define FACT_H

#include <gmp.h>
#include <stdio.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"

void print_results(const Params *P, const Results *R);

void contfract_factor(const Params *P, Results *R);

size_t choose_s_fb(const mpz_t N); 

#endif 
