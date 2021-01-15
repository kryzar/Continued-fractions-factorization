/* lp_var.c */

#include "lp_var.h"
#include "step_B.h"

int is_Qn_fact_lp_var(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows, const mpz_t Qn,
                      mpz_t lp, const mpz_t *factor_base, const size_t s_fb, 
                      const mpz_t pm_squared) {
    /* 
    Test (with trial division) if Qn completely factors over the factor
    base factor_base or almost completely factors (ie if after having 
    removed all prime divisors of Qn which belong to the factor base, 
    the remaining cofactor is less than pm_squared).
    During this process, when a prime factor of Qn is found, if this
    factor has an odd power in the factorization, its index in the
    factor base array is stored in the 'Qn_odd_pows' array. For each
    call of the 'is_Qn_fact_lp_var' function , the 'Qn_odd_pows' array 
    is filled (with dummy values) from the beginning. If Qn is factorisable
    or almost factorisable, this data is then used to create its exponent
    vector. The number of odd power is set in nb_Qn_odd_pows. 

    return: 1 if Qn is factorisable, -1 if Qn is almost completely 
            factorisable, 0 otherwise.
    param Qn_odd_pows: An array to store the indexes of the prime 
                       factors of Qn that have an odd valuation. It
                       is already allocated and initialized (size s_fb).
    param nb_Qn_odd_pows: A pointer to the number of indexes
                          added to the array (we start from the 
                          beginning). 
    param Qn: The Qn to be factored. 
    param lp: An auxiliary variable used to simplify Qn with the 
              primes of the factor base. If Qn factors almots completely, 
              it will store the large prime of the factorization.
    param factor_base: The factor base. 
    param s_fb: The size of the factor base.
    param pm_squared: The value pm^2 where pm is the largest prime of 
                      the factor base.
    */

    
    mp_bitcnt_t valuation; // Qn valuation for a prime of factor_base
    size_t i_fb ;     // Index of a prime of the factor_base array. 

    mpz_set(lp, Qn);  // We are going to simplify lp by the primes 
                      // of the factor base 
    *nb_Qn_odd_pows = 0; // To start from the beginning of the array
    i_fb = 0; 
     
    while (i_fb < s_fb &&  mpz_cmp_ui(lp, 1) ) {
        valuation = mpz_remove(lp, lp, factor_base[i_fb]); 
        if (1 == (valuation & 0x1) ) {
            Qn_odd_pows[*nb_Qn_odd_pows] = i_fb;
            (*nb_Qn_odd_pows) ++;  
        }
        i_fb ++; 
    }

    if ( 0 == mpz_cmp_ui(lp, 1) ) { 
        // If lp has been completely simplified
        return 1;  
    }else if (0 > mpz_cmp(lp, pm_squared) ) {
        // If the remaining cofactor is less than pm^2
        return -1; 
    }
    return 0;
}

AQp_lp *create_AQp_lp(const mpz_t Qn, const mpz_t Anm1, const mpz_t lp,
        exp_vect_data *D, size_t n) {
    /*
    This function allocates the memory for a node and initializes its
    member Qn, Anm1 and lp and exp_vect.

    param Qn: The Qn of a pair (Anm1, Qn) with Qn almost completely 
               factorisable.
    param Anm1: The Anm1 of a pair(Anm1, Qn) with Qn almost completely
                 factorisable.
    paramp lp: The large prime of the Qn's factorization.
    param D: A pointer to the structure containing the data to compute 
             the exponent_vector. 
    param n: The subscript of Qn needed to compute the exponent vector.
    */
    
    AQp_lp *node; 
      
    node = (AQp_lp *)malloc(sizeof(AQp_lp)); 
      
    mpz_init_set(node -> Qn, Qn); 
    mpz_init_set(node -> Anm1, Anm1); 
    mpz_init_set(node -> lp, lp); 
    init_exp_vect(1, node -> exp_vect, D, n); 
    node -> next = NULL; 

    return node; 
}

void delete_AQp_lp_list(AQp_lp **list) {
    /* 
    This function deletes the linked list. It deallocates its memory and
    sets its head pointer to NULL.
    
    param **list: A pointer to the head pointer of the list. 
    */
 
    AQp_lp *current; 
    AQp_lp *node;     // Node to be deleted

    current = *list; 
    while(current) {
        node = current; 
        current = current -> next; 
        mpz_clears(node -> Anm1, node -> Qn, node -> lp, node -> exp_vect, NULL);
        free(node); 
    }
    *list = NULL; 
}

void insert_or_elim_lp(AQp_lp **list, const mpz_t Qn, const mpz_t Anm1,
                    const mpz_t lp, exp_vect_data *D, size_t n, mpz_t *Qns,
                    mpz_t *Ans, mpz_t *exp_vects, size_t *nb_AQp, 
                    const mpz_t N, mpz_t A, mpz_t Q, mpz_t gcd, mpz_t exp_vect,
                    Results *R) {
    /* 
    Given a large prime lp, this fonction checks if lp is already a 
    parameter of a node of list, a sorted linked list. If it is, the AQp_lp
    with this large prime is used as the pivot of the Gaussian elimination 
    to eliminate the large prime and to add a new A-Q pair the Ans, Qns 
    and exp_vect arrays. If not, a node is created. After the elimination 
    process, if we have a new A-Q pair with Qn a square, it may be possible 
    to find a non trivial factor of N. If so, the factor is stored in R-> factor
    and R-> found is set to 1.
  
    param list: A pointer to the head pointer of the linked list.
    param Qn: The Qn of a pair (Anm1, Qn) with Qn almost completely factorisable.
    param Anm1: The Anm1 of a pair (Anm1, Qn) with Qn almost completely factorisable.
    param lp: The large prime of Qn.
    param D: The data needed to compute the exponent vector associated to Qn.
    param n: The subscript of Qn.
    param Qns: Is used if a new A-Q pair with Qn completely factorisable is 
               computed.
    param Ans: idem
    param exp_vects: idem
    param nb_AQp: idem
    param N: The number to be factored.
    param A: An auxiliary variable already initialized to perform one step 
             of the gaussian elimination.
    param Q: idem
    param gcd: An auxiliary variable already initialized used to see if we can
               find a non trivial factor of N if we have a new A-Q pair with Qn
               a square.
    param exp_vect: An auxiliary variable already initialized to perform one step 
                    of the gaussian elimination.
    param R: A pointer to the structure used to store, if we can,a non trivial
             factor of N.
    */
       
    AQp_lp *current; 

    current = *list;

    if (current == NULL) { 
        // If the list is empty, create a node and return it as the head
        // of the list.
        *list = create_AQp_lp(Qn, Anm1, lp, D, n);
        return;
    }
    else if (0 < mpz_cmp(current -> lp, lp)) { 
        // If the value of lp is smaller than the lp of the head node,
        // Create a node and insert it at the start of the list.
        AQp_lp *node = create_AQp_lp(Qn, Anm1, lp, D, n); 
        node -> next = current;
        *list = node; 
        return; 
    }

    while(current -> next) {
        if (0 < mpz_cmp(current -> next -> lp, lp)) { 
            if(mpz_cmp(current -> lp, lp)){
                // There is not yet in the linked list an AQp_lp with   
                // this value of lp : create a AQp_lp and insert it       
                AQp_lp *node = create_AQp_lp(Qn, Anm1, lp, D, n); 
                node -> next = current -> next; 
                current -> next = node; 
                return; 
            }else{
                // Do not add a nod. Use the Aqp_lp in the list which
                // has the same lp as the pivot of the gaussian elim.
                // Add the resulting A-Q pair to the Qns, Ans and 
                // exp_vects arrays.
                                       
                
                mpz_mul(Q, Qn, current  -> Qn);     // Multiply the 2 Qi.   
                mpz_mul(A, Anm1, current -> Anm1);  // Multiply the 2 Ai.
                mpz_mod(A, A, N);

                init_exp_vect(0, exp_vect, D, n); // Xor the 2 exponent vectors.
                mpz_xor(exp_vect, exp_vect, current -> exp_vect); 

                if (0 == mpz_cmp_ui(exp_vect, 0)) { 
                    // If Q is a square, ie if exp_vect == 0 
                     mpz_sqrt(gcd, Q);        // gcd <-- gcd(A - sqrt(Q), N)
                     mpz_sub(gcd, A, gcd);  
                     mpz_gcd(gcd, gcd, N);    
                    if (mpz_cmp_ui(gcd, 1) && mpz_cmp(gcd, N)) {
                        // We found a non trivial factor of N 
                        mpz_set(R->fact_found, gcd); 
                        R->found = 1; 
                        return;
                     }
                }else{
                    // The pair A-Q is a pair such that all the primes 
                    // that have an odd valuation in the factorization of Q
                    // belong to the factor base 
                  mpz_init_set(Qns[*nb_AQp], Q); 
                  mpz_init_set(Ans[*nb_AQp], A);
                  mpz_init_set(exp_vects[*nb_AQp], exp_vect); 
                  (*nb_AQp) ++;
                }
                return; 
            }
        }
        current = current -> next;
    }

    // Create a AQp_lp and insert in at the end of the list
    AQp_lp *node = create_AQp_lp(Qn, Anm1, lp, D, n); 
    current -> next = node; 
    return; 

}


void create_AQ_pairs_lp_var(const Params *P, Results *R, mpz_t *Ans, mpz_t *Qns, 
                            mpz_t *exp_vects, const mpz_t *factor_base,
                            AQp_lp **list) {
    /*
    This function computes the A-Q pairs, by expanding sqrt(kN) into a
    continued fraction and stores them in Ans and Qns. If Qn is completely
    factorisable with the primes of the factor base, it directly stores
    the A-Q pair in the Qns and Ans array and adds the exponent vector
    associated to Qn in the exp_vects array. If Qn is almost completely
    factorisable (see the large prime) variation, the auxiliary function
    insert_or_elim_lp is called to see if its large prime lp has already
    been encountered. If that is the case, the large prime is present in
    list and can be eliminated, performing one step of the gaussian
    elimination and the resulting A-Q pair is added in the Ans, Qns and 
    exp_vects array. If not, a struct AQp_lp is added to the sorted linked
    list list.
    If a Qn is a square, it may be possible to find a factor of N. In 
    this case, the factor is set in R->factor and R->found is set to 1. 

    param P: Pointer to the set of parameters for the problem (see step_A.h)
    param R: Pointer to the structure used to store the result (see step_A.h).
             (the structure isn't already initialized)
    param Ans: Array of size P->nb_want_AQp (already allocated but not 
               initialized) to store the An's.
    param Qns: Array of size P->nb_want_AQp (already allocated but not
               initialized) to store the Qn's.
    param exp_vects: Array of size P->nb_want_AQp (already allocated
                     but not initialized) to store the exponent vectors.
    param factor_base: The factor base of size P->s_fb. 
    param list: A linked list (NULL at the beginning). (see the 
                description of the struc AQp_lp in lp_var.h)
    */

    /***************
	* Declarations *
	***************/
 
    // For the auxiliary functions
    struct exp_vect_data D; 
    mpz_t  pm_squared; // To store the square of the largest prime of factor_base.
    mpz_t  lp;         // To store, if Qn has a large prime, that large prime.
    mpz_t  A;          // Auxiliary variable for insert_or_elim_lp
    mpz_t  Q;          // Auxiliary variable for insert_or_elim_lp
    mpz_t  gcd;        // Auxiliary variable for insert_or_elim_lp
    mpz_t  exp_vect;   // Auxiliary variable for insert_or_elim_lp
    int    r;          // To store the result of the 'is_Qn_fact_lp_var' and
                       // 'insert_or_elim_lp' functions

 
	D.Qn_odd_pows = (size_t *)malloc(P->s_fb * sizeof(size_t)); 
    D.reduced_fb_indexes = (size_t *)malloc(P->s_fb * sizeof(size_t));
    mpz_inits(pm_squared, lp, A, Q, gcd, exp_vect, NULL);

    mpz_mul(pm_squared, factor_base[P->s_fb - 1], factor_base[P->s_fb - 1]); 
    D.nb_reduced_fb_indexes = 0;
    
    // For the continued fraction expansion
    mpz_t  Anm1; // A_{n-1}
    mpz_t  An; 
    mpz_t  Qnm1; // Q_{n-1}
    mpz_t  Qn; 
    mpz_t  rnm1; // r_{n-1}
    mpz_t  rn; 
    mpz_t  qn; 
    mpz_t  Gn; 
    mpz_t  g; 
    mpz_t  temp; 
    mpz_t  AQtemp; // To store temporarily a An or Qn value
    size_t n;      // The subscript n in Qn
    size_t nb_AQp; // Number of found A-Q pairs found with Qn factorisable
                   // over the factor base. 
       
    mpz_inits(An, Qnm1, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL); 
    mpz_init_set_ui(Anm1, 1);        // A_{-1} <-- 1
    mpz_init_set_ui(Qn, 1);          // Q0 <-- 1

    mpz_mul_ui(Qnm1, P->N, P->k);      // Q_{-1} <-- kN
    mpz_sqrt(g, Qnm1);               // g = [sqrt(k*N)] 
    mpz_set(An, g);                  // A0 <-- g = [sqrt(k*N)]
    mpz_set(rnm1, g);                // r_{-1} <-- g
    mpz_set(qn, g);                  // q0 <-- g
    n      = 0;  
    nb_AQp = 0; 
 
    while ( n < P->n_lim && nb_AQp < P->nb_want_AQp) {
   
        /************
		* Expansion *
		************/

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
        mpz_mod(An, An, P->N);
        mpz_set(Anm1, AQtemp); 

        n++;
        
        // Is Qn factorisable ? 
        r = is_Qn_fact_lp_var(D.Qn_odd_pows, &(D.nb_Qn_odd_pows), Qn, lp,
                              factor_base, P->s_fb, pm_squared); 
   
        if (1 == r) { 
            // If Qn is completely factorisable
            if ( !(n & 0x1) && (0 == D.nb_Qn_odd_pows) ) {
                // If Qn is a square with n even: Anm1^2 = sqrt(Qn)^2 mod N.
                mpz_sqrt(temp, Qn);        // temp <-- sqrt(Qn)
                mpz_sub(temp, Anm1, temp); // temp <-- Anm1 - sqrt(Qn)
                mpz_gcd(temp, temp, P->N);  // temp <-- gcd(Anm1 - sqrt(Qn), N)
                // We may find a non trivial factor of N        
                if (mpz_cmp_ui(temp, 1) && mpz_cmp(temp, P->N)) { 
                    mpz_set(R->fact_found, temp);
                    R->found = 1; 

                    free(D.Qn_odd_pows); D.Qn_odd_pows = NULL; 
                    free(D.reduced_fb_indexes); D.reduced_fb_indexes = NULL; 
                    mpz_clears(pm_squared, lp, A, Q, gcd, exp_vect, Anm1, An,
                               Qnm1, Qn, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL); 
                    return; 
                }
            }else{
                // If the exponent vector associated to Qn is not zero
                mpz_init_set(Ans[nb_AQp], Anm1); // Store A_{n-1}
                mpz_init_set(Qns[nb_AQp], Qn);   // Store Qn
                init_exp_vect(1, exp_vects[nb_AQp], &D, n); 
                nb_AQp++; 
            }

        }else if (-1 == r){ 
            // If Qn is almost completely factorisable
            insert_or_elim_lp(list, Qn, Anm1, lp, &D, n, Qns, Ans, exp_vects,
                                &nb_AQp, P->N, A, Q, gcd, exp_vect, R);
            if (R->found){
                free(D.Qn_odd_pows); D.Qn_odd_pows = NULL; 
                free(D.reduced_fb_indexes); D.reduced_fb_indexes = NULL; 
                mpz_clears(pm_squared, lp, A, Q, gcd, exp_vect, Anm1, An, Qnm1, Qn, rnm1,
                           rn, qn, Gn, g, temp, AQtemp, NULL); 
                return; 
            }
        }
    }
    
    R->n_last = n; 
    R->nb_AQp = nb_AQp; 

	/*******
	* Free *
	*******/
    free(D.Qn_odd_pows); D.Qn_odd_pows = NULL; 
    free(D.reduced_fb_indexes); D.reduced_fb_indexes = NULL; 
    mpz_clears(pm_squared, lp, A, Q, gcd, exp_vect, Anm1, An, Qnm1, Qn, rnm1, rn,
               qn, Gn, g, temp, AQtemp, NULL);    
} 

