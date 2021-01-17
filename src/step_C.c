/* step_C.c */

#include <gmp.h>
#include "step_A.h"
#include "step_C.h"

void gauss_elimination(mpz_t *exp_vects, mpz_t *hist_vects, size_t *lin_rel_indexes,
	                 size_t *nb_lin_rel, const size_t nb_AQp) {
    /*
    Perform the gaussian elimination to determine if the exponents 
    vectors in exp_vects are linearly dependent. If so, the array 
    'lin_rel_indexes' stores the indexes of the lines of the matrix
    where a linear relation was found. The history vectors in hist_vects
    allow to keep track of the lines which were xored (they form the 
    right side matrix when performing the gaussian elimination by hand).

    param exp_vects: An array which contains the exponent vectors.
                     Each exponent vector represents a line of the
                     matrix. Calculations are performed directly on
                     the vectors.
    param hist_vects: An array which contains the history vectors.
                      At the beginning, they form the identity matrix.
                      Calculations are performed directly on the
                      vectors.
    param lin_rel_indexes: An nb_AQp-size array already allocated that
                           will contain the indexes of the history
                           vectors associated to an acceptable A-Q pairs
                           set.
    param nb_lin_rel: Pointer to the number of linear relations found.
    param nb_AQp: The size of the exp_vects, hist_vects and 
                  lin_rel_indexes arrays.
    */

    size_t msb_indexes[nb_AQp]; // msb_indexes[i] is the index (numbered
                                // from 1 because of the return of the 
                                // mpz_sizeinbase function) of the most
                                // significant bit of exp_vects[i]. If
                                // exp_vects[i] is 0, msb_indexes[i] = 0.
    size_t j; // Index (numbered from 1) of the column being processed.
    size_t pivot; // It will contain the smallest i such that the most
                  // significant bit of exp_vects[i] is in the jth column.

    *nb_lin_rel = 0; 

    // Init msb_indexes and j
    msb_indexes[0] = mpz_sizeinbase(exp_vects[0], 2);
    j = msb_indexes[0];
    
    for (size_t i = 1; i < nb_AQp; i++) {
        msb_indexes[i] = mpz_sizeinbase(exp_vects[i], 2);
        if (j < msb_indexes[i]) {
            j = msb_indexes[i];
        }
    }

    /**********************
    * Reduction procedure *
    **********************/

    while (j >= 1) {
        // Find the pivot if it exists
        pivot = 0;
        while (pivot < nb_AQp && msb_indexes[pivot] != j) {
            pivot ++;
        }

        if (pivot < nb_AQp) { // If we find a pivot
            for (size_t i = pivot + 1; i < nb_AQp; i++) {
                // Search for other i such that the most
                // significant bit of exp_vects[i] is in the
                // jth column.
                if (msb_indexes[i] == j) {
                    // Reduction
                    mpz_xor(hist_vects[i], hist_vects[i], hist_vects[pivot]);
                    mpz_xor(exp_vects[i], exp_vects[i], exp_vects[pivot]); 

                    if (0 == mpz_cmp_ui(exp_vects[i], 0)) {
                        // If we find a linear relation
                        lin_rel_indexes[*nb_lin_rel] = i;
                        (*nb_lin_rel) ++;
                        msb_indexes[i] = 0;
                    } else {
                        msb_indexes[i] = mpz_sizeinbase(exp_vects[i], 2);
                    }
                }
            }
        }
        j --;
    }
}


void calculate_A_Q(mpz_t A, const mpz_t *Ans, mpz_t Q, const mpz_t *Qns,
                   mpz_t hist_vect, const mpz_t N) {
    /*
    Find the A and the Q of the S-congruence A^2 = Q^2 mod N and put
    the results in A and Q.

    param A: Already initialized, used to store the value of A.
    param Ans: An array which contains the A_n's.
    param Q: Already initialized, used to store the value of Q.
    param Qns: An array which contains the Q_n's. 
    param hist_vect: An history vector which is associated to an S-Set.
                    (at the end, the vector is zero).
    param N: The integer to be factored.
   	*/

    // In the comments, let Q1, Q2, ... be the Q_i of the S-Set and
    // A1, A2, ... be the A_i of the S-Set.
    
    mp_bitcnt_t i;
    mpz_t R; 
    mpz_t X; 
    mpz_t Q_temp;

    mpz_inits(R, X, Q_temp, NULL); 
     
    mpz_set_ui(Q, 1);			 // Q <- 1
    i = mpz_scan1(hist_vect, 0);
    mpz_clrbit(hist_vect, i);
    mpz_set(A, Ans[i]);			 // A <- A1
    mpz_set(R, Qns[i]);			 // R <- Q1

    while  (0 != mpz_cmp_ui(hist_vect, 0)) {
        i = mpz_scan1(hist_vect, i+1);
        mpz_clrbit(hist_vect, i);

        mpz_mul(A, A, Ans[i]); // A <- A * Ai
        mpz_mod(A, A, N);

		mpz_set(Q_temp, Qns[i]); // X <- pgcd(R, Q_i)
        mpz_gcd(X, R, Q_temp);

        mpz_mul(Q, Q, X); // Q <- QX mod N
        mpz_mod(Q, Q, N);

        mpz_divexact(Q_temp, Q_temp, X); // R <- R / X * Q_i / X
        mpz_divexact(R, R, X);
        mpz_mul(R, R, Q_temp);
    }

    mpz_sqrt(X, R); // X <- sqrt(R)
    mpz_mul(Q, Q, X); // Q <- QX mod N 
    mpz_mod(Q, Q, N); 

    mpz_clears(R, X, Q_temp, NULL); 
}

void find_factor(Results *Res, const mpz_t *Ans, const mpz_t *Qns, 
                 mpz_t *exp_vects, mpz_t *hist_vects, const mpz_t N) {
    /*
    Find the S-Sets from the A-Q pairs previously calculated using the
    auxiliary function gaussian_elimination. After the reduction procedure,
    the 1's (plural) of a history vector whose index is in the array 
    lin_rel_indexes indicate the indexes of A-Q pairs in a S-Set.
    For each S-Set, find A and Q such that A^2 = Q^2 mod N using the
    auxiliary function calculate_A_Q. Compute the pgcd(A-Q, N) hoping
    to find a factor of N. If a factor is found, put it in R-> fact_found
    and set R-> found to 1.

    params Res: Pointer to the structure used to store the results.
    param Ans: The Ans computed by create_AQ_pairs 
    params Qns: The Qns computed by create_AQ_pairs
    param exp_vects: The exponent vectors computed by create_AQ_pairs
    param hist_vects: The historys vectors computed by init_hist_vects.
    param N: The integer to be factored
    */

    size_t *lin_rel_indexes; 
    size_t nb_lin_rel;
    mpz_t  A; 
    mpz_t  Q;   
    mpz_t  gcd;

    lin_rel_indexes = (size_t *)malloc(Res-> nb_AQp * sizeof(size_t)); 
    mpz_inits(A, Q, gcd, NULL);

    gauss_elimination(exp_vects, hist_vects, lin_rel_indexes, &nb_lin_rel,
					  Res-> nb_AQp);

    for (size_t i = 0; i < nb_lin_rel; i++) {
        calculate_A_Q(A, Ans, Q, Qns, hist_vects[lin_rel_indexes[i]], N);
        mpz_sub(gcd, A, Q);    // gcd <- pgcd(A - Q, N)
        mpz_gcd(gcd, gcd, N);  
        if (mpz_cmp_ui(gcd, 1) && mpz_cmp(gcd, N)) {
            mpz_set(Res-> fact_found, gcd);
            Res-> found = 1;

            free(lin_rel_indexes); lin_rel_indexes = NULL;
            mpz_clears(A, Q, gcd, NULL);
            return;
        }
    }

    free(lin_rel_indexes); lin_rel_indexes = NULL;
    mpz_clears(A, Q, gcd, NULL);

}
