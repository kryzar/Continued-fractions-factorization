/* main.c */

#include <time.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"
#include "fact.h"
#include "test.h"

int main(int argc, char **argv) {

    /*
    char *file_name = "graph_nb_bits_time";  
    nb_bits_vs_time(file_name, 80, 140, 5);  
    */

    Params P;
    Results R;

    init_Params_Results(&P, &R);
     
    mpz_set_str(P.N, "340282366920938463463374607431768211457", 10);
     
    P.k = choose_k(P.N, K_MAX); 
    P.n_lim = 3300000; 
   
    /*
    P.lp_var = 0;
    P.nb_want_AQp = 2060; 
    P.s_fb = 2700;
    */
    
    P.lp_var = 1; 
    P.s_fb = 450;
    P.nb_want_AQp = 450 + 16; 

    set_eas_params(&P, EAS_CUT, EAS_COEFF);
    contfract_factor(&P, &R); 
    print_results(&P, &R); 
    clear_Params_Results(&P, &R); 

    return 0; 
}
