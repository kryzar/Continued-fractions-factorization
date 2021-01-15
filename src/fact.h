/* fact.h */

#ifndef FACT
#define FACT

#include <gmp.h>
#include <stdio.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"
#include "lp_var.h"

void print_results(const Params *P, const Results *R);

void contfract_factor(const Params *P, Results *R); 

#endif 
