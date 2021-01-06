/*
	lp_var.c
*/

#include "lp_var.h"
#include "step_B.h"

int is_qn_fact_lp_var(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows, const mpz_t Qn,
                      mpz_t lp, const mpz_t *factor_base, const size_t s_fb, 
                      const mpz_t pm_squared){

      /* This function tests by trial division if Qn is completely 
       * factorisable with the primes of the factor base or almost
       * completely factorisable (ie if after having removed all prime
       * divisors of Qn which belong to the factor base, the remaining
       * cofactor is less than pm^2, where pm is the largest prime of 
       * the factor base).
       * During this process, when a prime factor of Qn is found, if this
       * factor has an odd power in the factorization, its index in the
       * factor base array is stored in the 'Qn_odd_pows' array. For each
       * call of the 'is_qn_fact_lp_var' function , the 'Qn_odd_pows' array 
       * is filled from the beginning. If Qn is factorisable or almost 
       * factorisable with the primes of the gactor base, these data
       * are then used to create its exponent vector. 
       */

      /* return: 1 if Qn is factorisable.
       *        -1 if Qn is almost completely factorisable.
       *         0 otherwise.
       * param Qn_odd_pows: An array to store the indexes of the prime 
       *                    factors of Qn that have an odd valuation. It
       *                    is already allocated and initialized (size s_fb).
       * param nb_Qn_odd_pows: A pointer to the number of indexes
       *                      added to the array (we start from the 
       *                      beginning). 
       * param Qn: The Qn to be factored. 
       * param lp: An auxiliary variable used to simplify Qn with the 
       *           primes of the factor base. If Qn is almost completely
       *           factorisable, it will store the large prime of the 
       *           factorization.
       * param factor_base: The factor base. 
       * param s_fb: The size of the factor base.
       * param pm_squared: The value pm^2 where pm is the largest prime 
       *                   of the factor base.
       */

      mp_bitcnt_t valuation; 
      size_t i_fb ;     // Index of a prime of the factor_base array. 

      mpz_set(lp, Qn);  // We are going to simplify lp by the
                        //primes of the factor base 
      *nb_Qn_odd_pows = 0; // To start from the beginning of the array
      i_fb = 0; 
     
      while (i_fb < s_fb &&  mpz_cmp_ui(lp, 1) ){
            valuation = mpz_remove(lp, lp, factor_base[i_fb]); 
            if (1 == (valuation & 0x1) ){
                  Qn_odd_pows[*nb_Qn_odd_pows] = i_fb; 
                  (*nb_Qn_odd_pows) ++;  
            }
            i_fb ++; 
      }
      if ( 0 == mpz_cmp_ui(lp, 1) ){ 
            // If lp has been completely simplified
            return 1;  
      }else if (0 > mpz_cmp(lp, pm_squared)){
            // If the remaining cofactor is less than pm^2
            return -1; 
      }
      return 0;
}

AQp_lp *create_AQp_lp(const mpz_t Qn, const mpz_t Anm1, const mpz_t lp, Data_exp_vect *D_exp_vect, size_t n){

      /* This function allocates the memory for a node and initializes 
       * its member Qn, Anm1 and lp and exp_vect. */

      /* param Qn: The Qn of a pair (Anm1, Qn) with Qn almost completely factorisable.
       * param Anm1: The Anm1 of a pair(Anm1, Qn) with Qn almost completely factorisable.
       * paramp lp: The large prime of the Qn's factorization.
       * param D_exp_vect: A pointer to the structure containing the data to compute the exponent_vector. 
       * param n: The subscript of Qn needed to compute the exponent vector.
       */

      AQp_lp *node; 
      
      node = (AQp_lp *)malloc(sizeof(AQp_lp)); 
      
      mpz_init_set(node -> Qn, Qn); 
      mpz_init_set(node -> Anm1, Anm1); 
      mpz_init_set(node -> lp, lp); 
      init_exp_vect(1, node -> exp_vect, D_exp_vect, n); 
      node -> next = NULL; 

      return node; 
}

void delete_AQp_lp_list(AQp_lp **list){

      /*  This function deletes the linked list. It deallocates its memory
       *  and sets its head pointer to NULL. */
       
      // param **list: a pointer to the head pointer of the list. 
       

      AQp_lp *current; 
      AQp_lp *node;     // Node to be deleted

      current = *list; 
      while(current){
            node = current; 
            current = current -> next; 
            mpz_clears(node -> Anm1, node -> Qn, node -> lp, node -> exp_vect, NULL);
            free(node); 
      }
      *list = NULL; 
}

AQp_lp *insert_or_eliminate_lp(AQp_lp *list, const mpz_t Qn, const mpz_t Anm1,
                              const mpz_t lp, Data_exp_vect *D_exp_vect, size_t n, mpz_t *Qns,
                              mpz_t *Ans, mpz_t *exp_vects, size_t *nb_AQp, 
                              const mpz_t N, mpz_t A, mpz_t Q, mpz_t gcd, mpz_t exp_vect){
      
      /* Given a large prime lp, this fonction checks if lp is already 
       * a parameter of a node of the sorted linked list 'list'. If it is, 
       * the AQp_lp with this large prime is used as the pivot of the Gaussian
       * elimination to eliminate the large prime and add a new A-Q pair to
       * the Ans, Qns and exp_vect arrays. If not, a node is created. */
       
      /* return: The head pointer of the modified linked list.
       * param list: The head pointer of the linked list.
       * param Qn: The Qn of a pair (Anm1, Qn) with Qn almost completely factorisable.
       * param Anm1: The Anm1 of a pair (Anm1, Qn) with Qn almost completely factorisable.
       * param lp: The large prime of Qn.
       * param D_exp_vect: The data needed to compute the exponent vector associated to Qn.
       * param n: The subscript of Qn.
       * param Qns: Is used if a new A-Q pair with Qn completely factorisable is computed.
       * param Ans: idem
       * param exp_vects: idem
       * param nb_AQp: idem
       * param N: The number to be factored.
       * param A: An auxiliary variable already initialized to perform one step 
       *          of the gaussian elimination.
       * param Q: idem
       * param gcd: An auxiliary variable already initialized used to see if we
       *            can find a non trivial factor of N if we have a new A-Q pair
       *            with Qn a square.
       * param exp_vect: An auxiliary variable already initialized to perform one step 
       *                 of the gaussian elimination.
       */

      AQp_lp *current; 

      current = list;

      if (current == NULL){                     // If the list is empty
            // Create a node and return it as the head of the list.
            return create_AQp_lp(Qn, Anm1, lp, D_exp_vect, n); 
      }
      else if (0 < mpz_cmp(current -> lp, lp)){ // If the value of lp is smaller 
                                                // than the lp of the head node:
            // Create a node and insert it at the start of the list.
            AQp_lp *node = create_AQp_lp(Qn, Anm1, lp, D_exp_vect, n); 
            node -> next = current; 
            return node; 
      }
      while(current -> next){
            if (0 < mpz_cmp(current -> next -> lp, lp)){ 
                  if(mpz_cmp(current -> lp, lp)){

                        /******************************************************
                         * There is not yet in the linked list an AQp_lp with  * 
                         * this value of lp : create a AQp_lp and insert it   *
                         ******************************************************/
                        AQp_lp *node = create_AQp_lp(Qn, Anm1, lp, D_exp_vect, n); 
                        node -> next = current -> next; 
                        current -> next = node; 
                        return list; 

                  }else{
                        /***********************************************************
                         * There is in the list an AQp_lp whith the same lp. (do   *
                         * not add a node). Instead, use its data to eliminate     *
                         * the large prime before adding the resulting A-Q pair to *
                         * the Qns, Ans and exp_vects arrays.                      *
                         **********************************************************/ 
                        
                        // One step of the gaussian elimination 
                        mpz_mul(Q, Qn, current  -> Qn);     // Multiply the 2 Qi.  
                        mpz_mod(Q, Q, N); 
                        mpz_mul(A, Anm1, current -> Anm1);  // Multiply the 2 Ai.
                        mpz_mod(A, A, N);

                        init_exp_vect(0, exp_vect, D_exp_vect, n); // Xor the 2 exponent vectors.
                        mpz_xor(exp_vect, exp_vect, current -> exp_vect); 

                        if (0 == mpz_cmp_ui(exp_vect, 0)){ // If Q is a square, ie if exp_vect == 0 
                              mpz_sqrt(gcd, Q);        // gcd <-- gcd(A - sqrt(Q), N)
                              mpz_sub(gcd, A, gcd);  
                              mpz_gcd(gcd, gcd, N);    
                              if (mpz_cmp_ui(gcd, 1) && mpz_cmp(gcd, N)){
                                     // A CHANGER : Que faire quand on trouve un facteur ? 
                                     gmp_printf("factor (found during the elimination of a large prime) : %Zd \n", gcd);   
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
                        return list; 
                  }
            }
            current = current -> next;
      }
      // Create a AQp_lp and insert in at the end of the list
      AQp_lp *node = create_AQp_lp(Qn, Anm1, lp, D_exp_vect, n); 
      current -> next = node; 
      return list; 

}

int create_AQ_pairs_lp_var(const Params P, mpz_t *Ans, mpz_t *Qns, size_t *nb_AQp,
                            mpz_t *exp_vects, const mpz_t *factor_base, 
                            const size_t s_fb, AQp_lp **list, mpz_t fact_found){
      
      /* This function computes, by expanding sqrt(kN) into a continued
       * fraction, the A-Q pairs.
       * ie a pair (A_{n-1}, Qn) such that A_{n-1}^2 = (-1)^n * Qn mod N.
       * If Qn is completely factorisable with the primes of the factor base,
       * it directly stores the A-Q pair in the Qns and Ans array and adds 
       * the exponent vector associated to Qn in the exp_vects array.
       * If Qn is almost completely factorisable (see the large prime) variation,
       * the auxiliary function 'insert_or_eliminate' is called to see if its 
       * large prime lp has already been encountered. If that is the case, the 
       * large prime is present in 'list' and can be eliminated, performing one
       * step of the gaussian elimination. We then add a pair A-Q in 
       * the Ans, Qns and exp_vects array. If not, a struct AQp_lp is added to 
       * the sorted linked list 'list'.
       */

      /* return: 1 if a non trivial factor was found.
       *         0 otherwise.
       * Param P: Set of parameters for the problem (see step_A.h)
       * Param Ans: Array of size P.nb_want_AQp (already allocated
       *             but not initialized) to store the An's.
       * Param Qns: Array of size P.nb_want_AQp (already allocated
       *             but not initialized) to store the Qn's.
       * Param nb_AQp: A pointer to the number of A-Q pairs found 
       *                with Qn factorisable with the primes of the
       *                factor base or almost factorisable.
       * Param exp_vects: Array of size P.nb_want_AQp (already allocated
       *                   but not initialized) to store the exponent vectors.
       * Param factor_base: The factor base.
       * Param s_fb: The size of the factor_base array.
       * Param list: A linked list (NULL at the beginning). (see the 
       *              description of the struc AQp_lp in lp_var.h)
       * Param fact_found: A mpz_t already initialized to store, if 
       *                   we can, a non trivial factor of N.
       */             
       
           
     /***************************************************************************
	* Declarations, allocations and initializations for the auxilary functions *
	***************************************************************************/

      struct Data_exp_vect D; 
      mpz_t pm_squared; // To store pm^2 where pm is the largest prime of the factor base.
      mpz_t lp;         // To store, if Qn has a large prime, that large prime.
      mpz_t A;          // Auxiliary variable for insert_or_eliminate_lp
      mpz_t Q;          // Auxiliary variable for insert_or_eliminate_lp
      mpz_t gcd;        // Auxiliary variable for insert_or_eliminate_lp
      mpz_t exp_vect;   // Auxiliary variable for insert_or_eliminate_lp
      int r;            // To store the result of the 'is_qn_fact_lp_var' function

	D.Qn_odd_pows = (size_t *)malloc(s_fb * sizeof(size_t)); 
      D.reduced_fb_indexes = (size_t *)malloc(s_fb * sizeof(size_t));
      mpz_inits(pm_squared, lp, A, Q, gcd, exp_vect, NULL);

      mpz_mul(pm_squared, factor_base[s_fb - 1], factor_base[s_fb - 1]); 
      D.nb_reduced_fb_indexes = 0;

      /****************************************************************************
	*  Declarations, initializ. and assignments for the cont. frac. expansion   *
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
      size_t n; 
      
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

            r = is_qn_fact_lp_var(D.Qn_odd_pows, &(D.nb_Qn_odd_pows), Qn, lp, factor_base, s_fb, pm_squared); 

            if (1 == r){ // If Qn is completely factorisable
                  if ( !(n & 0x1) && (0 == D.nb_Qn_odd_pows) ){
                        // If Qn is a square with n even: Anm1^2 = sqrt(Qn)^2 mod N.
                        mpz_sqrt(temp, Qn);        // temp <-- sqrt(Qn)
                        mpz_sub(temp, Anm1, temp); // temp <-- Anm1 - sqrt(Qn)
                        mpz_gcd(temp, temp, P.N);  // temp <-- gcd(Anm1 - sqrt(Qn), N)
                        if (mpz_cmp_ui(temp, 1) && mpz_cmp(temp, P.N)){ // If we find a non trivial factor
                              mpz_set(fact_found, temp); 
                              return 1; 
                        }
                  }else{
                        // If the exponent vector associated to Qn is not zero
                        mpz_init_set(Ans[*nb_AQp], Anm1); // Store A_{n-1}
                        mpz_init_set(Qns[*nb_AQp], Qn);   // Store Qn
                        init_exp_vect(1, exp_vects[*nb_AQp], &D, n); 
                        (*nb_AQp)++; 
                  }
            }else if (-1 == r){ // If Qn is almost completely factorisable
                  *list = insert_or_eliminate_lp(*list, Qn, Anm1, lp, &D, n, Qns,
                              Ans, exp_vects, nb_AQp, P.N, A, Q, gcd, exp_vect); 
            }
      }

	/*****************
	* Free and clear *
	******************/
      free(D.Qn_odd_pows); D.Qn_odd_pows = NULL; 
      free(D.reduced_fb_indexes); D.reduced_fb_indexes = NULL; 
      mpz_clears(pm_squared, lp, A, Q, gcd, exp_vect, Anm1, An, Qnm1, Qn, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL);

      return 0; 
 
}
