/* step_B.c */

#include "step_B.h"
    
int is_Qn_factorisable(const Params *P, size_t *Qn_odd_pows,
                       size_t *nb_Qn_odd_pows, const mpz_t Qn, mpz_t Qn_divided,
                       const mpz_t *factor_base, const mpz_t pm_squared) {

    /*
    Test (with trial division) if Qn is factor_base-smooth. If the large
    prime variation is used it also tests if Qn is almost factor_base-smooth
    (if after having removed all prime divisors of Qn which belong to
    factor_base, the remaining cofactor is less than pm_squared). If the 
    early abort strategy is used, give up Qn if after P-> eas_cut divisions
    the unfactored prtion of Qn exceeds P-> eas_bound_div.
    During this process, when a prime factor of Qn is found and if this
    factor has an odd power in the factorization, its index in factor_base 
    is stored in Qn_odd_pows. For each call of the function is_Qn_factorisable,
    the array Qn_odd_pows is filled (with dummy values) from the beginning. 
    If Qn is factorisable, the dummy values become usefull: this data is then
    used to create its exponent vector. The number of odd powers is set in 
    nb_Qn_odd_pows.

    return: 1 if Qn is factor_base-smooth, -1 if the large prime variation
            is used and Qn is almost factor_base-smooth, 0 otherwise.
    param P: Set of parameters from the problem used to know the size of
             factor_base and if the large prime variation should be used.
    param Qn_odd_pows: An array to store the indexes of the prime 
                       factors of Qn that have an odd valuation. It
                       is already allocated and initialized (size s_fb).
    param nb_Qn_odd_pows: A pointer to the number of indexes
                          added to the array (we start from the
                          beginning).
    param Qn: The Qn to be factored.
    param Qn_divided: An auxiliary variable used to simplify Qn with
                      the primes of factor_base. If Qn is almost
                      factor_base-smooth it will store at the end 
                      the large prime of Qn's the factorization.
    param factor_base: The factor base.
    param pm_squared: The value pm^2 where pm is the largest prime of 
                      factor_base (only used with the large prime 
                      variation)
    */
    
    mp_bitcnt_t valuation;   // Qn valuation for a prime of factor_base
    mpz_set(Qn_divided, Qn); // We are going to divide Qn_divided by the
                             // primes of the factor base
    *nb_Qn_odd_pows = 0;     // To start from the beginning of the array

    size_t i = 0;
    while (i < P-> s_fb && mpz_cmp_ui(Qn_divided, 1) ) {
        // Early abort strategy 
        if (P-> eas && i == P-> eas_cut ) {
            if ( mpz_cmp(Qn_divided, P->eas_bound_div) > 0 ) {
                // If Qn_divided > P-> eas_bound_div give up Qn
                return 0; 
            }
        }
        // Find the valuation of factor_base[i] and simplify Qn_divided
        valuation = mpz_remove(Qn_divided, Qn_divided, factor_base[i]); 
        if (1 == (valuation & 0x1) ) {
            // If the valuation is odd
            Qn_odd_pows[*nb_Qn_odd_pows] = i; 
            (*nb_Qn_odd_pows) ++; 
        }
        i++; 
    }

    if ( 0 == mpz_cmp_ui(Qn_divided, 1) ) {
        // If Qn_divided has been completely simplified
        return 1; 
    } else if (P-> lp_var && 0 > mpz_cmp(Qn_divided, pm_squared) ) {
        // If the large prime variation is used and if the remaining 
        // cofactor is less than pm^2 
        return -1; 
    }
    return 0; 
} 

void init_hist_vects(mpz_t *hist_vects, size_t nb_AQp) {
    /*
    Prepare the hist_vects array for the gaussian elimination: put
    nb_AQp history vectors in the array (all the vectors needed).
    The initial value (set here) of a history vector v is 00..010...0
    where the 1 is in the column i (numbered from right to left from 0)
    such that hist_vect[i] = v.
    
    param hist_vects: The array of history vectors. Is is already 
                      allocated (size P.nb_want_AQp) but isn't
                      initialized.
    param nb_AQp: The number of AQ pairs found with Qn factor_base 
                  smooth.(or almost smooth, with possibly a squared
                  large prime in its decomposition if the large prime
                  variation is used) 
    */
    for (size_t i = 0; i < nb_AQp; i++) {
        mpz_init(hist_vects[i]);
        mpz_setbit(hist_vects[i], i);
    }
}


void init_exp_vect(int init, mpz_t exp_vect, const size_t *Qn_odd_pows, 
                   size_t nb_Qn_odd_pows, size_t n) {
    /*
    Initialize exp_vect if init == 1 and compute its value with the
    data stored in the array Qn_odd_odd_pows and the index n. If n
    is odd, set the parity bit of exp_vect to 1. The indexes stored
    in Qn_odd_pows indicate the bits to set to 1.

    param init: 1 if exp_vect needs to be initialized, 0 otherwise.
    param exp_vect: The exponent vector to compute.
    param Qn_odd_pows: Array which stores the indexes of the prime 
                       factors of Qn that have an odd valuation.
    param nb_Qn_odd_pows: The number of such prime factors of Qn.
    param n: The index n which is important for the "parity bit".
    */

    if (init) {
        mpz_init(exp_vect);
    } else {
        mpz_set_ui(exp_vect, 0);
    }

    if (n & 0x1) {
        // Set the "parity bit" 
        mpz_setbit(exp_vect, 0); // Set the least significant bit to 1
    }

    for (size_t i = 0; i < nb_Qn_odd_pows; i++) {
        // Set the other bits
        mpz_setbit(exp_vect, Qn_odd_pows[i] + 1); 
    }
} 

void create_AQ_pairs(const Params *P, Results *R, mpz_t *Ans, mpz_t *Qns, 
                     mpz_t *exp_vects, const mpz_t *factor_base, 
                     AQp_lp **list) {
    /*
    Compute the A-Q pairs, by expanding sqrt(kN) into a continued fraction.
    If Qn is factor_base-smooth, store directly the A-Q pair in the Ans, 
    Qns arrays and add in exp_vects the exponent vector associated to Qn.
    If Qn is almost factor_base-smooth and if the large prime variation
    is used, call the auxiliary function insert_or_elim_lp to see if its
    large prime has already been encountered. If that is the case, the 
    large prime is present in list and is eliminated, performing one step
    of the gaussian elimination and the resulting A-Q pair is added in the
    Ans, Qns and exp_vects array. If not, a struct AQp_lp is added to the
    sorted linked list list.
    Set the number of A-Q pairs found in R-> nb_AQp and the total number
    of A-Q pairs computed (factor_base-smooth or not)in R-> n_last.
    If a Qn is a square, it may be possible to find a factor of N. In this
    case, set the factor in R-> factor and set R-> found to 1.

    param P: Pointer to the set of parameters for the problem (see step_A.h).
    param R: Pointer to the structure used to store the result (see step_A.h).
             The structure is already initialized.
    param Ans: Array of size P-> nb_want_AQp (already allocated but not 
               initialized) to store the An's.
    param Qns: Array of size P-> nb_want_AQp (already allocated but not
               initialized) to store the Qn's.
    param exp_vects: Array of size P-> nb_want_AQp (already allocated
                     but not initialized) to store the exponent vectors.
    param factor_base: The factor base of size P-> s_fb. 
    param list: A linked list (NULL at the beginning). (see the 
                description of the structure in lp_var.h).
    */
    
    /********************************************************
    * Declarations for is_Qn_factorisable and init_exp_vect *
    ********************************************************/
    
    size_t *Qn_odd_pows;
    size_t  nb_Qn_odd_pows;     
    mpz_t   pm_squared; // To store the square of the largest prime of 
                        // factor_base.
    mpz_t   Qn_divided; // Used to simplify Qn with the primes of factor_base.
                        // If Qn is almost factor_base-smooth it will store 
                        // its large prime.
    
    Qn_odd_pows = (size_t *)malloc(P-> s_fb * sizeof(size_t));
    mpz_inits(pm_squared, Qn_divided, NULL); 

    mpz_mul(pm_squared, factor_base[P-> s_fb - 1], factor_base[P-> s_fb - 1]);

    /*************************************
    * Declarations for insert_or_elim_lp *
    **************************************/

    mpz_t  A;          
    mpz_t  Q;          
    mpz_t  gcd;        
    mpz_t  exp_vect; 
                
    mpz_inits(A, Q, gcd, exp_vect, NULL);
 
    /***********************************************
     * Declarations for the cont. fract. expansion *
     **********************************************/
    mpz_t  An;
    mpz_t  Anm1;   // A_{n-1}
    mpz_t  Gn;
    mpz_t  Qn;
    mpz_t  Qnm1;   // Q_{n-1}
    mpz_t  g;
    mpz_t  qn;
    mpz_t  rn;
    mpz_t  rnm1;   // r_{n-1}
    mpz_t  temp;
    mpz_t  AQtemp; // To store temporarily an An or Qn value
    size_t n;	   // The subscript n in Qn
    size_t nb_AQp; // Number of A-Q pairs such that Qn has all
                   // its prime factor with an odd valuation in
                   // factor_base.
    int    r;      // To store the result of is_Qn_factorisable
    
    mpz_inits(An, Qnm1, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL);
    mpz_init_set_ui(Anm1, 1);       // A_{-1} <-- 1
    mpz_init_set_ui(Qn, 1);         // Q0 <-- 1

    mpz_mul_ui(Qnm1, P-> N, P-> k); // Q_{-1} <-- kN
    mpz_sqrt(g, Qnm1);              // g = [sqrt(k*N)] 
    mpz_set(An, g);                 // A0 <-- g = [sqrt(k*N)]
    mpz_set(rnm1, g);               // r_{-1} <-- g
    mpz_set(qn, g);                 // q0 <-- g
    n	    = 0;
    nb_AQp  = 0;

    /*************
     * Expansion *
     ************/

    while (n < P-> n_lim && nb_AQp < P-> nb_want_AQp) {

        // Perform one step of the expansion
 
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
        mpz_mod(An, An, P-> N);
        mpz_set(Anm1, AQtemp);

        n++;
        
        // Is Qn factorisable ?
        r = is_Qn_factorisable(P, Qn_odd_pows, &nb_Qn_odd_pows, Qn, 
                               Qn_divided, factor_base, pm_squared);
        if (1 == r) {
            // If Qn is factor_base-smooth
            if ( !(n & 0x1) && (0 == nb_Qn_odd_pows) ) {
                // If Qn is a square with n even: Anm1^2 = sqrt(Qn)^2 mod N
                // we may find a non trivial factor of N
                mpz_sqrt(temp, Qn);        // temp <-- sqrt(Qn)
                mpz_sub(temp, Anm1, temp); // temp <-- Anm1 - sqrt(Qn)
                mpz_gcd(temp, temp, P-> N);  // temp <-- gcd(Anm1 - sqrt(Qn), N)
                if (mpz_cmp_ui(temp, 1) && mpz_cmp(temp, P-> N)) { 
                    // We find a factor
                    mpz_set(R-> fact_found, temp);
                    R-> found = 1; 
                    free(Qn_odd_pows); Qn_odd_pows = NULL; 
                    mpz_clears(pm_squared, Qn_divided, A, Q, gcd, exp_vect, Anm1, An,
                               Qnm1, Qn, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL); 
                    return; 
                }
            } else {
                // If the exponent associated to Qn is not zero
                // Add a A-Q pair
                mpz_init_set(Ans[nb_AQp], Anm1); // Store A_{n-1}
                mpz_init_set(Qns[nb_AQp], Qn);   // Store Qn
                init_exp_vect(1, exp_vects[nb_AQp], Qn_odd_pows,
                              nb_Qn_odd_pows, n); 
                nb_AQp++; 
            }
        } else if (-1 == r) { 
            // If Qn is almost B-smooth and the large prime variation is used
            init_exp_vect(0, exp_vect, Qn_odd_pows, nb_Qn_odd_pows, n); 
            insert_or_elim_lp(list, Qn, Anm1, Qn_divided, exp_vect, Qns, Ans,
                              exp_vects, &nb_AQp, P-> N, A, Q, gcd, R); 
            if (R-> found) {
                free(Qn_odd_pows); Qn_odd_pows = NULL;
                mpz_clears(pm_squared, Qn_divided, A, Q, gcd, exp_vect, Anm1, An, Qnm1, Qn, rnm1,
                           rn, qn, Gn, g, temp, AQtemp, NULL);
                return;
            }
        }
    }

    R-> n_last = n;
    R-> nb_AQp = nb_AQp;

    /*******
    * Free *
    *******/

    free(Qn_odd_pows); Qn_odd_pows = NULL;
    mpz_clears(pm_squared, Qn_divided, A, Q, gcd, exp_vect, Anm1, An, Qnm1, Qn, rnm1, rn,
               qn, Gn, g, temp, AQtemp, NULL);
}
