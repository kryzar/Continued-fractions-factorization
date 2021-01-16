/* step_B.h */

#ifndef STEP_B_H
#define STEP_B_H

#include <gmp.h>
#include "step_A.h"

typedef struct exp_vect_data{
	/*
	This structure gathers all the data needed to compute the exponent
	vector of a given Qn, except the subscript n.
	*/

	size_t* reduced_fb_indexes;   /* The elements of this array indicate
									 the index of the primes of the
									 factor base that have an odd
									 valuation in at least one of the
									 Qn's factorization (with Qn
									 completely factored or almost
									 completely factored if the large
									 prime variation is used). */ 
	size_t nb_reduced_fb_indexes; /* The number of indexes added to the 
                                     reduced_fb_indexes array. */
	size_t* Qn_odd_pows;          /* Store the indexes of the prime
									 factors of Qn that have an odd
									 valuation. */
	size_t nb_Qn_odd_pows;        /* The number of such prime factors of
									 Qn. */
     
} exp_vect_data;

typedef struct AQp_lp{
    /*
    The goal of this structure is to create a linked list whose
    nodes store the data of an A-Q pairs (Anm1, Qn), where Qn 
    is almost completely factorisable with the primes of the factor
    base. The remaining cofactor is called lp (for large prime). 
    This linked list has two properties. It is stored according
    to the member lp. Only one AQp_lp exists for a given value of
    lp. That AQp_lp pair is used as the pivot of the Gaussian 
    elimination to eliminate the large prime lp if another pair
    with the same large prime cofactor is encounterded. 
     */
    mpz_t Qn; 
    mpz_t Anm1; 
    mpz_t lp;               // The large prime.  
    mpz_t exp_vect;         // The exponent vector associated to Qn
                            // (the large prime lp is missing).
    struct AQp_lp *next;    // The next node of the linked list. 
}AQp_lp; 
 
int is_Qn_factorisable(const Params *P, size_t *Qn_odd_pows, 
                       size_t *nb_Qn_odd_pows, const mpz_t Qn, mpz_t Qn_divided,
                       const mpz_t *factor_base, const mpz_t pm_squared); 

void init_hist_vects(mpz_t *hist_vects, const size_t nb_AQp);

void init_exp_vect(const int init, mpz_t exp_vect, exp_vect_data *D,
				   const size_t n); 

void create_AQ_pairs(const Params *P, Results *R, mpz_t *Ans, mpz_t *Qns, 
                     mpz_t *exp_vects, const mpz_t *factor_base, 
                     AQp_lp **list);

AQp_lp *create_AQp_lp(const mpz_t Qn, const mpz_t Anm1, const mpz_t lp, 
                      exp_vect_data *D, size_t n); 

void insert_or_elim_lp(AQp_lp **list, const mpz_t Qn, const mpz_t Anm1,
                              const mpz_t lp, exp_vect_data *D, size_t n,
                              mpz_t *Qns, mpz_t *Ans, mpz_t *exp_vects, 
                              size_t *nb_AQp, const mpz_t N, mpz_t A, mpz_t Q,
                              mpz_t gcd, mpz_t exp_vect, Results *R);

void delete_AQp_lp_list(AQp_lp **list);

#endif
