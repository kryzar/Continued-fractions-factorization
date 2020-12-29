#include "step_C.h"
#include <gmp.h>


void gauss_elimination(mpz_t *exp_vects, mpz_t *hist_vects, size_t *lin_rel_indexes,
                       size_t *nb_lin_rel, const size_t nb_AQp){



      /* This function performs the gaussian elimination to determine 
       * whether the set of exponent vectors is linearly dependent.
       * If so, the array 'lin_rel_indexes' will store the indexes of 
       * the lines of the matrix where a linear relation was found. The
       * history vectors keep track of the lines which were xored.
       */

      /* param exp_vects: An array which contains the exponent vectors.
       *                  Each exponent vector represents a line of the
       *                  matrix. Calculus are directly performed on 
       *                  the vectors.
       * param hist_vects: An array which contains the history vectors.
       *                   At the beginning, they form the identity
       *                   matrix. Calculus are directly performed on
       *                   the vectors. 
       * param lin_rel_indexes: An array already allocated of size 
       *                        'nb_AQp'. It will contain the indexes 
       *                        of the history vectors associated to an
       *                        S-set. 
       * param nb_lin_rel: a pointer to the number of linear relations
       *                     found.
       * param nb_AQp: The size of the exp_vects, hist_vects and 
       *               lin_rel_indexes arrays. 
       */

      size_t msb_indexes[nb_AQp]; /* msb_indexes[i] is the index (numbered
                                     from 1) of the most significant  
                                     bit of exp_vects[i]. If exp_vects[i] 
                                     is zero, msb_indexes[i] = 0.
                                     */
      size_t i; 
      size_t j; 
      size_t pivot; /* It will contain the smallest i such that the 
                       most significant bit of exp_vects[i] is in the
                       jth column. */

      *nb_lin_rel = 0; 

      /* Initialization of msb_indexes
         and j, the index (numbered from 1) of the leftmost column */
      msb_indexes[0] = mpz_sizeinbase(exp_vects[0], 2); 
      j = msb_indexes[0]; 
      for (i = 1; i < nb_AQp; i++){
            msb_indexes[i] = mpz_sizeinbase(exp_vects[i], 2); 
            if (j < msb_indexes[i]){
                  j = msb_indexes[i]; 
            }
      }

      // Reduction procedure 
      
      while ( j >= 1 ){

            // Find the pivot if it exists 
            pivot = 0; 
            while (pivot < nb_AQp && msb_indexes[pivot] != j){
                  pivot ++; 
            }

            if (pivot < nb_AQp){ // If we find a pivot
                  for (i = pivot + 1; i < nb_AQp; i++){     
                        /* Search for other i such that the most significant
                        bit of exp_vects[i] is in th jth column*/ 
                        if (msb_indexes[i] == j){ 
                              // Reduction 
                              mpz_xor(hist_vects[i], hist_vects[i], hist_vects[pivot]); 
                              mpz_xor(exp_vects[i], exp_vects[i], exp_vects[pivot]); 
                        
                              if (0 == mpz_cmp_ui(exp_vects[i], 0)){ 
                                    // If we find a linear relation
                                    lin_rel_indexes[*nb_lin_rel] = i; 
                                    (*nb_lin_rel) ++;
                                    msb_indexes[i] = 0;
                              }else{
                                    msb_indexes[i] = mpz_sizeinbase(exp_vects[i], 2); 
                              }

                        }
                  }
            }
            j --; 
      }
}


void find_A_Q(mpz_t A, const mpz_t *Ans, mpz_t Q, const mpz_t *Qns, mpz_t hist_vect,
              const mpz_t N, mpz_t R, mpz_t X, mpz_t Q_temp){

      /* This function finds the A and the Q of the S-congruence A^2 = Q^2 mod N.
       
       * param A: Is already initialized. It is used to store the value of A.
       * param Ans: An array which contains the A_n's.
       * param Q: Is already initialized. It is used to store the value of Q.
       * param Qns: An array which contains the Q_n's. 
       * param hist_vect: An history vector which is associated to an S-Set.
       *                  (at the end, the vector is zero).
       * param N: The integer to be factored.
       * param R: An auxiliary variable already initialized for the computation of Q. 
       * param X: An auxiliary variable already initialized for the computation of Q.
       * param Q_temp: An auxiliary variable already initialized for the computation of Q. 
       */

       /* In the comments, let Q1, Q2, ... be the Q_i of the S-Set 
        *  and A1, A2, ... be the A_i of the S-Set */

      mp_bitcnt_t i;

      mpz_set_ui(Q, 1); // Q <-- 1

      i = mpz_scan1(hist_vect, 0);
      mpz_clrbit(hist_vect, i);

      mpz_set(A, Ans[i]); // A <-- A1
      mpz_set(R, Qns[i]); // R <-- Q1

      while  ( 0 != mpz_cmp_ui(hist_vect, 0) ){
            i = mpz_scan1(hist_vect, i+1);
            mpz_clrbit(hist_vect, i);

            mpz_mul(A, A, Ans[i]); // A <-- A * Ai 
            mpz_mod(A, A, N);
  
            mpz_set(Q_temp, Qns[i]); // X <-- pgcd(R, Q_i)
            mpz_gcd(X, R, Q_temp);

            mpz_mul(Q, Q, X); // Q <-- QX mod N 
            mpz_mod(Q, Q, N); 
 
            mpz_divexact(Q_temp, Q_temp, X); // R <-- R / X * Q_temp / X
            mpz_divexact(R, R, X); 
            mpz_mul(R, R, Q_temp);
      }
      mpz_sqrt(X, R); // X <-- sqrt(R)
      mpz_mul(Q, Q, X); // Q <-- QX mod N 
      mpz_mod(Q, Q, N); 
}



void find_A(mpz_t A, const mpz_t *Ans, mpz_t hist_vect, const mpz_t N){
      // @Antoine, comme les calculs de A et de Q suivent exactement la 
      // même structure de boucle, boucle qui nécessite de trouver 
      // le bit à 1 de hist_vect et de le modifier, il vaut mieux les faire
      // dans une même fonction. Je mets ici la fonction pour le calcul de 
      // A seul et en dessous la fonction pour le calcul de Q seul. Fonctions
      // que j'ai fusionnées dans find_A_Q. 
      mp_bitcnt_t i;

      i = mpz_scan1(hist_vect, 0);
      mpz_clrbit(hist_vect, i);
      mpz_set(A, Ans[i]); 

      while  ( 0 != mpz_cmp_ui(hist_vect, 0) ){
            i = mpz_scan1(hist_vect, i+1);
            mpz_clrbit(hist_vect, i); 
            mpz_mul(A, A, Ans[i]); 
            mpz_mod(A, A, N);
      }
} 


void find_Q(mpz_t Q, const mpz_t *Qns, mpz_t hist_vect, const mpz_t N,
            mpz_t R, mpz_t X, mpz_t Q_temp){
      // à ne pas utiliser 
      mp_bitcnt_t i; 

      mpz_set_ui(Q, 1);
      
      i = mpz_scan1(hist_vect, 0); // R <-- Q_1 
      mpz_clrbit(hist_vect, i); 
      mpz_set(R, Qns[i]);
      
      while  ( 0 != mpz_cmp_ui(hist_vect, 0) ){
            i = mpz_scan1(hist_vect, i+1); // X <-- pgcd(R, Q_i)
            mpz_clrbit(hist_vect, i); // with Q1, Q2, ... the Q_i of the S-Set

            mpz_set(Q_temp, Qns[i]); 
            mpz_gcd(X, R, Q_temp);

            mpz_mul(Q, Q, X); // Q <-- QX mod N 
            mpz_mod(Q, Q, N); 
 
            mpz_divexact(Q_temp, Q_temp, X); // R <-- R / X * Q_temp / X
            mpz_divexact(R, R, X); 
            mpz_mul(R, R, Q_temp);
      }
      mpz_sqrt(X, R); // X <-- sqrt(R)
      mpz_mul(Q, Q, X); // Q <-- QX mod N 
      mpz_mod(Q, Q, N); 
}
