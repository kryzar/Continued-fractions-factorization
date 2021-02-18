/* main.c */

#include <gmp.h>
#include <time.h>
#include <string.h>
#include "init_algo.h"
#include "step_1.h"
#include "step_2.h"
#include "fact.h"
#include "test.h"

int main(int argc, char **argv) {

    mpz_t       N; 
    mp_bitcnt_t nb_bits; 

    if (argc != 2) {
        printf("Entrez un argument en ligne de commande (cf README.txt). \n"); 
    } else if ( 0 == strcmp(argv[1], "F7_lp_eas") ) {
        fact_F7(1,1); 
    }  else if ( 0 == strcmp(argv[1], "F7_lp") ) {
        fact_F7(1,0); 
    } else if ( 0 == strcmp(argv[1], "F7_eas") ) {
        fact_F7(0,1); 
    } else if ( 0 == strcmp(argv[1], "F7") ) {
        fact_F7(0,0); 
    } else {
        mpz_init_set_str(N, argv[1], 10); 
        if ( mpz_cmp_ui(N, 150) > 0) {
            fact_N(N); 
        } else {
            nb_bits = mpz_get_ui(N); 
            fact_rand_N(nb_bits);
        }
        mpz_clear(N); 
    }

    return 0; 
}
