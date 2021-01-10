/* main.c */

#include <stdio.h>
#include <time.h>
#include "step_A.h"
#include "step_B.h"
#include "step_C.h"
#include "lp_var.h"

int main(int argc, char **argv) {

      Params P; 
      mpz_t *Ans;
      mpz_t *Qns; 
      mpz_t *exp_vects; 
      mpz_t *hist_vects; 
      mpz_t *factor_base;
      mpz_t fact_found; 
      size_t s_fb; 
      size_t nb_AQp; 
      size_t i;
      int found; 
      AQp_lp *list_AQp_lp = NULL; // For the large prime variation 

       /************************
       *  Just for the test    *
       *************************/

      mpz_init_set_str(P.N, "340282366920938463463374607431768211457", 10); 
      P.n_lim =  1330000;
      P.k = 257;

      //P.nb_want_AQp = 2060; 
      //s_fb = 2700; 
      
      P.nb_want_AQp = 650; 
      s_fb = 700; 
   
       /************************
       *      Allocations      *
       *************************/

      Ans = malloc_mpz_array(P.nb_want_AQp); 
      Qns = malloc_mpz_array(P.nb_want_AQp);
      exp_vects = malloc_mpz_array(P.nb_want_AQp); 
      hist_vects = malloc_mpz_array(P.nb_want_AQp);
      factor_base = malloc_mpz_array(s_fb);
      mpz_init(fact_found); 

       /************************
       *  Looking for a factor  *
       *************************/
      init_factor_base(factor_base, s_fb, P.N, P.k);
      found = create_AQ_pairs_lp_var(P, Ans, Qns, &nb_AQp, exp_vects, factor_base, 
                                    s_fb, &list_AQp_lp, fact_found);                            // for the large prime variation

      //found = create_AQ_pairs(P, Ans, Qns, &nb_AQp, exp_vects, factor_base, s_fb, fact_found); // without the large prime variation 
      
      if (! found){
            init_hist_vects(hist_vects, nb_AQp);
            found = find_factor(Ans, Qns, exp_vects, hist_vects, nb_AQp, P.N, fact_found); 
      }
                  
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
      if (found){
            gmp_printf("factor: %Zd \n", fact_found); 
      }else{
            printf("factor not found \n"); 
      }

       /************************
       *          Free         *
       *************************/
      free_mpz_array(Ans, nb_AQp); 
      free_mpz_array(Qns, nb_AQp); 
      free_mpz_array(exp_vects, nb_AQp); 
      free_mpz_array(hist_vects, nb_AQp); 
      free_mpz_array(factor_base, s_fb); 
      delete_AQp_lp_list(&list_AQp_lp); // for the large prime variation
      mpz_clear(fact_found); 

      return 0; 
}
