/* lp_var.c */

#include "lp_var.h"

AQp_lp *create_AQp_lp(const mpz_t Qn, const mpz_t Anm1, const mpz_t lp, 
                      const mpz_t exp_vect) {
    /*
    Allocate the memory for a node and initialize its member Qn, Anm1,
    lp and exp_vect.

    return: The node created. 
    param Qn: The Qn of a pair (Anm1, Qn) with Qn almost factor_base
              smooth.
    param Anm1: The Anm1 of the A-Q pair. 
    paramp lp: The large prime of Qn's factorization.
    param exp_vect: The exponent vector associated to Qn.
    */

    AQp_lp *node = (AQp_lp *)malloc(sizeof(AQp_lp));

    mpz_init_set(node-> Qn, Qn); 
    mpz_init_set(node-> Anm1, Anm1); 
    mpz_init_set(node-> lp, lp);
    mpz_init_set(node-> exp_vect, exp_vect); 
    node-> next = NULL; 

    return node; 
} 

void delete_AQp_lp_list(AQp_lp **list) {
    /* 
    Delete the linked list: deallocate its memory and set its head
    pointer to NULL.
    
    param **list: A pointer to the head pointer of the list. 
    */
 
    AQp_lp *current; 
    AQp_lp *node;     // Node to be deleted

    current = *list; 
    while (current) {
        node = current; 
        current = current-> next; 
        mpz_clears(node-> Anm1, node-> Qn, node-> lp, node-> exp_vect, NULL);
        free(node); 
    }
    *list = NULL; 
}
    
void insert_or_elim_lp(AQp_lp **list, const mpz_t Qn, const mpz_t Anm1, 
                       const mpz_t lp, mpz_t exp_vect, mpz_t *Qns,
                       mpz_t *Ans, mpz_t *exp_vects, size_t *nb_AQp,
                       const mpz_t N, mpz_t A, mpz_t Q, mpz_t gcd,
                       Results *R) {
    /*
    Given a large prime lp, check if lp is already a parameter of a node
    of list (a sorted linked list). If it is, use the AQp_lp with this
    large prime as the pivot of the gaussian elimination to eliminate lp
    and to add a new A-Q pair in the Ans, Qns and exp_vects arrays. If
    not, create a node. 
    After the elimination process, if we have a new A-Q pair with Qn a
    square, it may be possible to find a non trivial factor of N. If so, 
    store the factor in R-> factor and set R-> found to 1. 

    param list: A pointer to the head pointer of the linked list. 
    param Qn: The Qn of a pair (Anm1, Qn) with Qn almost factor_base smooth.
    param Anm1: The Anm1 of the A-Q pair. 
    paramp lp: The large prime of Qn's factorization.
    param exp_vect: The exponent vector associated to Qn.
    param Qns: Used (if lp is already in the list) to store the new A-Q pair
               after having eliminated the large prime.
    param Ans: idem 
    param exp_vects: idem
    param nb_AQp: Pointer to number of A-Q pairs in Ans, Qns and exp_vects
                  (to be incremented if we add a new A-Q pair.)
    param N: The number to be factored.
    param A: Auxiliary variable already initialized to perform one step 
             of the gaussian elimination.
    param Q: Auxiliary variable already initialized to perform one step 
             of the gaussian elimination.
    param gcd: An auxiliary variable already initialized used to see if
               we can find a factor of N if we have a new A-Q pair with
               Qn a square.
    param R: A pointer to the structure used to store, if we can, a non
             trivial factor of N.    
    */

    AQp_lp *current;
    AQp_lp *node; // If a node must be created 

    current = *list;

    if (current == NULL) {
        // If the list is empty, create a node and make it the head of the list
        *list = create_AQp_lp(Qn, Anm1, lp, exp_vect);
        return; 
    } else if (0 < mpz_cmp(current-> lp, lp) ) {
        // If the value of lp is smaller than the lp of the head node
        // Create a node and insert it at the start of the list
        node = create_AQp_lp(Qn, Anm1, lp, exp_vect);
        node-> next = current; 
        *list = node; 
        return; 
    }

    while (current-> next) {
        if(0 < mpz_cmp(current-> next-> lp, lp)) {
            // We reach the position where the node may be inserted.
            if(mpz_cmp(current-> lp, lp)) {
                // There is not yet in the linked list an AQp_lp with 
                // this value of lp : create a node and insert it
                node = create_AQp_lp(Qn, Anm1, lp, exp_vect); 
                node-> next= current-> next;
                current-> next = node;
                return; 
            } else {
                // Do not add a nod. Use the Aqp_lp in the list which
                // has the same lp as a pivot of the gaussian elim.
                // Add the resulting A-Q pair to the Qns, Ans and 
                // exp_vects arrays.
    
                // Mulitply the 2 Qi and the 2 Ai
                mpz_mul(Q, Qn, current-> Qn);      
                mpz_mul(A, Anm1, current-> Anm1); 
                mpz_mod(A, A, N);
                // Xor the 2 exponent vectors.
                mpz_xor(exp_vect, exp_vect, current-> exp_vect);
                if (0 == mpz_cmp_ui(exp_vect, 0)) { 
                    // If Q is a square, ie if exp_vect == 0 
                    mpz_sqrt(gcd, Q);        // gcd <-- gcd(A - sqrt(Q), N)
                    mpz_sub(gcd, A, gcd);
                    mpz_gcd(gcd, gcd, N);    
                    if (mpz_cmp_ui(gcd, 1) && mpz_cmp(gcd, N)) {
                        // We found a non trivial factor of N 
                        mpz_set(R-> fact_found, gcd); 
                        R-> found = 1;
                    }
                } else {
                    // Add a new A-Q pair
                    mpz_init_set(Qns[*nb_AQp], Q); 
                    mpz_init_set(Ans[*nb_AQp], A);
                    mpz_init_set(exp_vects[*nb_AQp], exp_vect); 
                    (*nb_AQp) ++;
                }
                return; 
            }
        }
        current = current-> next;  
    }
    // Create a node and insert in at the end of the list
    node = create_AQp_lp(Qn, Anm1, lp, exp_vect);
    current-> next = node; 
}
