/* main.c */

#include <time.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"
#include "fact.h"

int main(int argc, char **argv) {

    Params P;
    Results R;

    mpz_init_set_str(P.N, "340282366920938463463374607431768211457", 10);
    P.k = choose_k(P.N, K_MAX); 
    P.n_lim = 3300000; 
   
    /*
    P.lp_var = 0;
    P.nb_want_AQp = 2060; 
    P.s_fb = 2700;
    */
    
    P.lp_var = 1; 
    P.s_fb = choose_s_fb(P.N);
    P.nb_want_AQp = P.s_fb + 16; 

    init_results(&R); 
    contfract_factor(&P, &R); 
    print_results(&P, &R); 
    clear_Params_Results(&P, &R); 

    return 0; 
}
