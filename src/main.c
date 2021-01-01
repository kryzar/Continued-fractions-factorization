/* main.c */

#include <stdio.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"

int main(int argc, char **argv) {

      Params P; 
      mpz_t *Ans;
      mpz_t *Qns; 
      mpz_t *exp_vects; 
      mpz_t *hist_vects; 
      mpz_t *factor_base;
      size_t s_fb; 
      size_t nb_AQp; 
      size_t i;

       /************************
       *  Just for the test    *
       *************************/
 
      mpz_init_set_str(P.N, "340282366920938463463374607431768211457", 10); 
      P.n_lim =  1400000;
      P.k = 257;
      P.nb_want_AQp = 2060; 
      s_fb = 2700; 
      
       /************************
       *      Allocations      *
       *************************/

      Ans = malloc_mpz_array(P.nb_want_AQp); 
      Qns = malloc_mpz_array(P.nb_want_AQp);
      exp_vects = malloc_mpz_array(P.nb_want_AQp); 
      hist_vects = malloc_mpz_array(P.nb_want_AQp);
      factor_base = malloc_mpz_array(s_fb);

       /************************
       *  Looking for a factor  *
       *************************/
      init_factor_base(factor_base, s_fb, P.N, P.k);
      create_AQ_pairs(P, Ans, Qns, &nb_AQp, exp_vects, hist_vects, factor_base, s_fb);
      find_factor(Ans, Qns, exp_vects, hist_vects, nb_AQp, P.N); 

       /************************
       *  Just for the test    *
       *************************/ 
      gmp_printf("N: %Zd \n", P.N); 
      printf("k: %u \n", P.k); 
      printf("s_fb : %lu \n", s_fb); 
      printf("n_lim: %lu \n", P.n_lim); 
      printf("nb_want_AQp: %lu \n", P.nb_want_AQp); 
      printf("nb_AQp: %lu \n", nb_AQp); 
      gmp_printf("dernier premier de la factor base: %Zd \n", factor_base[s_fb -1]); 


       /************************
       *          Free         *
       *************************/
      free_mpz_array(Ans, nb_AQp); 
      free_mpz_array(Qns, nb_AQp); 
      free_mpz_array(exp_vects, nb_AQp); 
      free_mpz_array(hist_vects, nb_AQp); 

      return 0; 
}
