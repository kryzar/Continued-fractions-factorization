/* main.c */

#include <stdio.h>
#include <time.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"
#include "lp_var.h"
#include "fact.h"

int main(int argc, char **argv) {

	Params P;
    Results R;

    mpz_init_set_str(P.N, "340282366920938463463374607431768211457", 10);
    P.k = 257; 
    P.n_lim = 1330000; 

    /*
    P.lp_var = 0;
    P.nb_want_AQp = 2060; 
    P.s_fb = 2700; 
    */

    P.lp_var = 1; 
    P.nb_want_AQp = 650; 
    P.s_fb = 700;

    init_results(&R); 
    contfract_factor(&P, &R); 
    print_results(&P, &R); 
    clear_Params_Results(&P, &R); 

    return 0; 
}
