/*
	step_B.c
*/

#include <gmp.h>
#include "step_A.h"
#include "step_B.h"


int is_qn_factorisable(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows, const mpz_t Qn,
                       mpz_t Q_temp, const mpz_t *factor_base, const size_t s_fb){

      /* This function tests by trial division if Qn is completely 
       * factorisable with the primes of the factor base. During this
       * process, when a prime factor of Qn is found, if this factor 
       * has an odd power in the factorization, its index in the factor
       * base array is stored in the 'Qn_odd_pows' array. For each call
       * of the 'is_qn_factorisable' function , the 'Qn_odd_pows' array 
       * is filled from the beginning. If Qn is factorisable, these data
       * are then used to create its exponent vector. 
       */

      /* return: 1 if Qn is factorisable, 0 otherwise. 
       * param Qn_odd_pows: An array to store the indexes of the prime 
       *                    factors of Qn that have an odd valuation. It
       *                    is already allocated and initialized (size s_fb).
       * param nb_Qn_odd_pows: The number of indexes
       *                      added to the array (we start from the 
       *                      beginning). 
       * param Qn: The Qn to be factored. 
       * param Q_temp: An auxiliary variable. 
       * param factor_base: The factor base. 
       * param s_fb: The size of the factor base.
       */

      mp_bitcnt_t valuation; 
      size_t i_fb ; // Index of a prime of the factor_base array. 

      mpz_set(Q_temp, Qn); /* We are going to simplify Q_temp by the
                            primes of the factor base */
      *nb_Qn_odd_pows = 0; // To start from the beginning of the array
      i_fb = 0; 
     
      while (i_fb < s_fb &&  mpz_cmp_ui(Q_temp, 1) ){
            valuation = mpz_remove(Q_temp, Q_temp, factor_base[i_fb]); 
            if (1 == (valuation & 0x1) ){
                  Qn_odd_pows[*nb_Qn_odd_pows] = i_fb; 
                  (*nb_Qn_odd_pows) ++;  
            }
            i_fb ++; 
      }
      if ( 0 == mpz_cmp_ui(Q_temp, 1) ){ 
            // If Q_tem has been completely simplified
            return 1; 
      } 
      return 0; 
}



void add_hist_vect(mpz_t *hist_vects, const size_t nb_AQp){

      /* This function adds in hist_vects the (nb_AQp + 1) th 
       * history vector. ie 2^(nb_AQp) written in base 2. 
       */

      /* param hist_vects: The array of history vectors. Is is already 
       *                    allocated (size P.nb_want_AQp) but isn't 
       *                    initialized.
       * param nb_AQp: The index which indicates where to put the (nb_AQp + 1)the
       *                history vector.
       */

      mpz_init(hist_vects[nb_AQp]); 
      mpz_setbit(hist_vects[nb_AQp], nb_AQp); 
}


void add_exp_vect(mpz_t *exp_vects, size_t *reduced_fb_indexes,
                  size_t *nb_reduced_fb_indexes, const size_t *Qn_odd_pows,
                  const size_t nb_Qn_odd_pows, const size_t nb_AQp, const size_t n){
      
      /* This function adds in exp_vects[nb_AQp] the exponent vector
       * associated to Q_n. 
      
       * param exp_vects: The array of exponent vectors. It is already
       *                  allocated (size P.nb_want_AQp) but isn't 
       *                  initialized.
       *                   exp_vects[0], ..., exp_vects[nb_AQp - 1] have
       *                   already been computed. 
       * param reduced_fb_indexes: An array already allocated, which is 
       *                            the same size as the factor base. 
       *                            It won't be completely filled. Its
       *                            elements indicate the index of the 
       *                            primes of the factor base that have 
       *                            an odd valuation in at least one of
       *                            the Qn's factorization (with Qn 
       *                            completely factored). 
       * param nb_reduced_fb_indexes: The number of indexes added to the
       *                              reduced_fb_indexes array.
       * param Qn_odd_pows:  It stores the indexes of the prime 
       *                    factors of Qn that have an odd valuation.
       * param nb_Qn_odd_pows: The number of such prime factors of Qn.
       * param nb_AQp: The number of A-Q pairs with Qn factorisable
       *                already found. 
       * param n: The subscript of Qn, which is important for the "parity bit". 
       */

      size_t i; 
      size_t Qn_odd_pow; 
      size_t bit_index; 
      
      mpz_init(exp_vects[nb_AQp]); 
      
      // Set the "parity bit" which indicates the sign of (-1)^n in the
      // expression A_{n-1}^2 = (-1)^n Qn [mod n]
      if (n & 0x1){
            mpz_setbit(exp_vects[nb_AQp], 0); // Set the least significant bit to 1
      }

      for (i = 0; i < nb_Qn_odd_pows; i++){
            Qn_odd_pow = Qn_odd_pows[i]; 
            // Check if Qn_odd_pow is in reduced_fb_indexes
            bit_index = 0; 
            while (bit_index < *nb_reduced_fb_indexes && Qn_odd_pow != reduced_fb_indexes[bit_index] ){
                  bit_index ++; 
            }
            // Set the bit of the exponent vector to 1 in the column (bit_index)
            // where Qn_odd_pow has been found or add a column and put a 1 in it
            mpz_setbit(exp_vects[nb_AQp], bit_index + 1); 
            // If Qn_odd_pow isn't in the reduced_fb_indexes array, we add it
            if (bit_index == *nb_reduced_fb_indexes){
                  reduced_fb_indexes[*nb_reduced_fb_indexes] = Qn_odd_pow; 
                  (*nb_reduced_fb_indexes) ++; 
            }
      }
}


void create_AQ_pairs(const Params P, mpz_t *Ans, mpz_t *Qns, size_t *nb_AQp,
                     mpz_t *exp_vects, mpz_t *hist_vects, 
                     const mpz_t *factor_base, const unsigned int s_fb){


      /* This function computes, by expanding sqrt(kN) into a continued
       * fraction, the A-Q pairs  
       * ie pairs (A_{n-1}, Q_n) such that A_{n-1}^2 = (-1)^n * Q_n mod N.
       * It only stores a pair if Qn is factorisable with the primes of
       * the factor base. In this case, it also adds the corresponding
       * exponent and history vector of Qn.
       */

      /* :Param P: Set of parameters for the problem (see step_A.h).
       * :Param Ans: Array of size P.nb_want_AQp (already allocated
       *             but not initialized) to store the An's.
       * :Param Qns: Array of size P.nb_want_AQp (already allocated
       *             but not initialized) to store the Qn's.
       * :Param nb_AQp: Number of A-Q pairs found with Qn factorisable
       *                with the primes of the factor base.
       * :Param exp_vects: Array of size P.nb_want_AQp (already allocated
       *                   but not initialized) to store the exponent vectors.
       * :Param hist_vects: Array of size P.nb_want_AQp (already allocated
       *                    but not initialized) to store the history vectors.
       * :Param factor_base: The factor base.
       * :Param s_fb: The size of the factor_base array.
       */

     /***************************************************************************
	* Declarations, allocations and initializations for the auxilary functions *
	***************************************************************************/
	
	size_t *Qn_odd_pows; 
      size_t *reduced_fb_indexes; 
      size_t  nb_Qn_odd_pows;   
      size_t nb_reduced_fb_indexes; 

      Qn_odd_pows = (size_t *)malloc(s_fb * sizeof(size_t)); 
      reduced_fb_indexes = (size_t *)malloc(s_fb * sizeof(size_t));
      nb_reduced_fb_indexes = 0;

     /****************************************************************************
	* Declarations, initializ. and assignments for the cont. frac. expansion  *
	****************************************************************************/
      mpz_t Anm1; // A_{n-1}
      mpz_t An; 
      mpz_t Qnm1; // Q_{n-1}
      mpz_t Qn; 
      mpz_t rnm1; // r_{n-1}
      mpz_t rn; 
      mpz_t qn; 
      mpz_t Gn; 
      mpz_t g; 
      mpz_t temp; 
      mpz_t AQtemp; // To store temporarily a An or Qn value
      size_t n; // The subscript n in Qn 

      mpz_inits(An, Qnm1, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL); 
      mpz_init_set_ui(Anm1, 1);        // A_{-1} <-- 1
      mpz_init_set_ui(Qn, 1);          // Q0 <-- 1

      mpz_mul_ui(Qnm1, P.N, P.k);      // Q_{-1} <-- kN
      mpz_sqrt(g, Qnm1);               // g = [sqrt(k*N)] 
      mpz_set(An, g);                  // A0 <-- g = [sqrt(k*N)]
      mpz_set(rnm1, g);                // r_{-1} <-- g
      mpz_set(qn, g);                  // q0 <-- g 
      n = 0; 
      *nb_AQp = 0; 

      /**********************************************************************/

      while ( n < P.n_lim && *nb_AQp < P.nb_want_AQp ){

            /**********************************************
            *                 Expand sqrt(kn)             *
            **********************************************/
 
            // Q_n = Q_{n-2} + q_{n-1} (r_{n-1} - r_{n-2})
            mpz_set(AQtemp, Qn); 
            mpz_sub(temp, rn, rnm1); 
            mpz_set(Qn, Qnm1);
            mpz_addmul(Qn, qn, temp);
            mpz_set(Qnm1, AQtemp); 

            // G_{n} = 2*g - r_{n-1}
            mpz_mul_ui(Gn, g, 2);
            mpz_sub(Gn, Gn, rn);
            
            // q_n = [G_n / Q_n]    r_n = G_n - q_n * Q_n  
            mpz_set(rnm1, rn); 
            mpz_fdiv_qr(qn, rn, Gn, Qn);
                      
            // An = q_n * A{n-1} + A{n-2} mod N
            mpz_set(AQtemp, An);
            mpz_mul(An, qn, An); 
            mpz_add(An, An, Anm1);
            mpz_mod(An, An, P.N);
            mpz_set(Anm1, AQtemp); 

            n++;
 
            /**********************************************
            * is (Anm1, Qn) a pair with Qn factorisable ? *
            **********************************************/

            if (is_qn_factorisable(Qn_odd_pows, &nb_Qn_odd_pows, Qn, temp, factor_base, s_fb)){ 
                  mpz_init_set(Ans[*nb_AQp], Anm1); // Store A_{n-1}
                  mpz_init_set(Qns[*nb_AQp], Qn);   // Store Qn 
                  add_hist_vect(hist_vects, *nb_AQp); 
                  add_exp_vect(exp_vects, reduced_fb_indexes, &nb_reduced_fb_indexes, Qn_odd_pows, nb_Qn_odd_pows, *nb_AQp, n); 
                  (*nb_AQp)++; 
            }
      }

	/*****************
	* Free and clear *
	******************/
      free(Qn_odd_pows); Qn_odd_pows = NULL; 
      free(reduced_fb_indexes); reduced_fb_indexes = NULL; 
      mpz_clears(Anm1, An, Qnm1, Qn, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL); 

}
