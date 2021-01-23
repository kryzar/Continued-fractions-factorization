/* fact.h */

#ifndef FACT_H
#define FACT_H

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include <time.h> 
#include "init_algo.h"
#include "step_1.h"
#include "step_2.h"

void print_results(const Params *P, const Results *R);
void contfract_factor(const Params *P, Results *R);

#endif 
