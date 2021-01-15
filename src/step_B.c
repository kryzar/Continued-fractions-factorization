/* step_B.c */

#include <gmp.h>
#include "step_A.h"
#include "step_B.h"


int is_Qn_factorisable(size_t *Qn_odd_pows, size_t *nb_Qn_odd_pows,
					   const mpz_t Qn, mpz_t Qn_divided, const mpz_t *factor_base,
					   const size_t s_fb) {
	/*
	Test (with trial division) if Qn completely factors over the factor
	base factor_base. During this process, when a prime factor of Qn is
	found and if this factor has an odd power in the factorization, its
	index in the factor base array is stored in Qn_odd_pows. For each
	call of the function is_Qn_factorisable, the array Qn_odd_pows is
	filled (with dummy values) from the beginning. If Qn is
	factorisable, the dummy values are replaced in Qn_odd_pows by real
	values and this data is then used to create its exponent vector. The
	number of odd powers is set in nb_Qn_odd_pows.

	return: 1 if Qn is factorisable, 0 otherwise. 
	param Qn_odd_pows: An array to store the indexes of the prime 
	                   factors of Qn that have an odd valuation. It
	                   is already allocated and initialized (size s_fb).
	param nb_Qn_odd_pows: A pointer to the number of indexes
	                      added to the array (we start from the 
	                      beginning). 
	param Qn: The Qn to be factored. 
	param Qn_divided: An auxiliary variable. 
	param factor_base: The factor base. 
	param s_fb: The size of the factor base.
	*/

	mp_bitcnt_t valuation; // Qn valuation for a prime of factor_base

	mpz_set(Qn_divided, Qn); // We are going to divide Qn_divided by the
                             // primes of the factor base
	*nb_Qn_odd_pows = 0; // To start from the beginning of the array
     
	size_t i = 0;
	while (i < s_fb &&  mpz_cmp_ui(Qn_divided, 1) ) {
		valuation = mpz_remove(Qn_divided, Qn_divided, factor_base[i]); 
		if (1 == (valuation & 0x1) ) {
			Qn_odd_pows[*nb_Qn_odd_pows] = i; 
			(*nb_Qn_odd_pows) ++;  
		}
		i ++; 
	}

	if ( 0 == mpz_cmp_ui(Qn_divided, 1) ) { 
		// If Q_divided has been completely simplified.
		return 1; 
	} 
	return 0; 
}

void init_hist_vects(mpz_t *hist_vects, const size_t nb_AQp) {
	/*
	Prepare the hist_vects array for the gaussian elimination :it puts
    nb_AQp history vectors in the array (all the vectors needed). 
    The initial value (set here) of a history vector v is 00..010...0
    where the 1 is in the column i (numbered from right to left from 0)
    such that hist_vect[i] = v. 

	param hist_vects: The array of history vectors. Is is already 
	                  allocated (size P.nb_want_AQp) but isn't 
	                  initialized.
	param nb_AQp: The number of AQ pairs found with Qn completely
	              factorisable with the primes of the factor base.
	*/

	for (size_t i = 0; i < nb_AQp; i++) {
		mpz_init(hist_vects[i]); 
		mpz_setbit(hist_vects[i], i); 
	}
}

void init_exp_vect(const int init, mpz_t exp_vect, exp_vect_data *D,
				   const size_t n) {
	/*
	This function initializes exp_vect if init = 1 and computes its
	value with the data stored in the struct exp_vect_data and the
	subscript n.

	param init: 1 if exp_vect needs to be initialized, 0 otherwise.
	param exp_vect: The exponent vector to compute (not initialized).
	param D: A pointer to the structure containing the data required
	         to compute the exponent vector (see step_B.h).
	         P_exp -> reduced_fb_indexes is an array (already
	         allocated) of same size as the factor base. It won't be
			 completely filled.
	param n: The subscript n which is important for the "parity bit".
	*/

	size_t i; 
	size_t Qn_odd_pow; 
	size_t bit_index; 

	if (init) {
		mpz_init(exp_vect); 
	} else {
		mpz_set_ui(exp_vect, 0); 
	}
		
	// Set the "parity bit" which indicates the sign of (-1)^n in the
	// expression A_{n-1}^2 = (-1)^n Qn [mod n]
	if (n & 0x1) {
		mpz_setbit(exp_vect, 0); // Set the least significant bit to 1
	}

	for (i = 0; i < D -> nb_Qn_odd_pows; i++) {
		Qn_odd_pow = D -> Qn_odd_pows[i]; 
		 
		// Check if Qn_odd_pow is in reduced_fb_indexes
		bit_index = 0; 
		while (bit_index < D -> nb_reduced_fb_indexes
			   && Qn_odd_pow != D -> reduced_fb_indexes[bit_index]) {
			bit_index ++; 
		}

		// Set the bit of the exponent vector to 1 in the column
		// (bit_index) where Qn_odd_pow has been found or add a column
		// and put a 1 in it
		mpz_setbit(exp_vect, bit_index + 1); 
		// If Qn_odd_pow isn't in the reduced_fb_indexes array, add it
		if (bit_index == D -> nb_reduced_fb_indexes) {
			D -> reduced_fb_indexes[D -> nb_reduced_fb_indexes] = Qn_odd_pow; 
			(D -> nb_reduced_fb_indexes) ++; 
		}
	}
} 

int create_AQ_pairs(const Params P, mpz_t *Ans, mpz_t *Qns, size_t *nb_AQp,
					mpz_t *exp_vects, const mpz_t *factor_base,
					const size_t s_fb, mpz_t fact_found) {
	/*
	This function computes the A-Q pairs, by expanding sqrt(kN) into a
	continued fraction, and stores them in Ans and Qns. It only stores a
	pair if Qn is completely factorisable with the primes of the factor
	base. In this case, it also computes Qn's exponent vector and adds
	it to exp_vects. The number of A-Q pairs found is set in nb_AQp. If
    a Qn is a square, it may be possible to find a factor of N.
    Return 1 is a non trivial factor of N was found, 0 otherwise.

	return: 1 if a non trival factor of Qn was found, 0 otherwise.
	param P: Set of parameters for the problem (see step_A.h).
	param Ans: Array of size P.nb_want_AQp (already allocated but not
			   initialized) to store the An's.
	param Qns: Array of size P.nb_want_AQp (already allocated but not 
			   initialized) to store the Qn's.
	param nb_AQp: A pointer to the number of found A-Q pairs found with
				  Qn factorisable of the factor base.
	param exp_vects: Array of size P.nb_want_AQp (already allocated but
					 not initialized) to store the exponent vectors.
	param factor_base: The factor base.
	param s_fb: The size of the factor_base array.
	param fact_found: A mpz_t already initialized to store, if we can,
					  a non trivial factor of N.
	*/

	/***************
	* Declarations *
	***************/
	
	// For the auxilary functions
	struct exp_vect_data D; 
	D.Qn_odd_pows			= (size_t *)malloc(s_fb * sizeof(size_t)); 
	D.reduced_fb_indexes	= (size_t *)malloc(s_fb * sizeof(size_t));
	D.nb_reduced_fb_indexes = 0;

	// For the continued fraction expansion
	mpz_t  AQtemp; // To store temporarily an An or Qn value
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
	size_t n;	   // The subscript n in Qn 

	mpz_inits(An, Qnm1, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL); 
	mpz_init_set_ui(Anm1, 1); // A_{-1} <-- 1
	mpz_init_set_ui(Qn, 1);	  // Q0 <-- 1

	mpz_mul_ui(Qnm1, P.N, P.k); // Q_{-1} <-- kN
	mpz_sqrt(g, Qnm1);          // g = [sqrt(k*N)] 
	mpz_set(An, g);             // A0 <-- g = [sqrt(k*N)]
	mpz_set(rnm1, g);           // r_{-1} <-- g
	mpz_set(qn, g);             // q0 <-- g 
	n		= 0; 
    *nb_AQp = 0; 

	while (n < P.n_lim && *nb_AQp < P.nb_want_AQp) {

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
		mpz_mod(An, An, P.N);
		mpz_set(Anm1, AQtemp); 

		n++;

		// Is Qn factorisable?
		if (is_Qn_factorisable(D.Qn_odd_pows, &(D.nb_Qn_odd_pows), Qn, temp,
							   factor_base, s_fb)) {
			if (!(n & 0x1) && (0 == D.nb_Qn_odd_pows)) {
				// If Qn is a square with n even: Anm1^2 = sqrt(Qn)^2 mod N.
				mpz_sqrt(temp, Qn);		   // temp <- sqrt(Qn)
				mpz_sub(temp, Anm1, temp); // temp <- Anm1 - sqrt(Qn)
				mpz_gcd(temp, temp, P.N);  // temp <- gcd(Anm1 - sqrt(Qn), N)
				// We may find a non trivial factor of N
				if (mpz_cmp_ui(temp, 1) && mpz_cmp(temp, P.N)) {
					mpz_set(fact_found, temp);

                    free(D.Qn_odd_pows); D.Qn_odd_pows = NULL; 
                    free(D.reduced_fb_indexes); D.reduced_fb_indexes = NULL; 
                    mpz_clears(Anm1, An, Qnm1, Qn, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL);

					return 1; 
				}
			} else {
				// If the exponent vector associated to Qn is not zero
				mpz_init_set(Ans[*nb_AQp], Anm1); // Store A_{n-1}
				mpz_init_set(Qns[*nb_AQp], Qn);   // Store Qn
				init_exp_vect(1, exp_vects[*nb_AQp], &D, n); 
				(*nb_AQp)++; 
			}
		}
	}

	/*******
	* Free *
	*******/
	
	free(D.Qn_odd_pows); D.Qn_odd_pows = NULL; 
	free(D.reduced_fb_indexes); D.reduced_fb_indexes = NULL; 
	mpz_clears(Anm1, An, Qnm1, Qn, rnm1, rn, qn, Gn, g, temp, AQtemp, NULL);

	return 0;

}
