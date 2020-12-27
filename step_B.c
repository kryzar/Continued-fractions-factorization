/*
	step_B.c
*/

#include <gmp.h>

#include <step_A.h>
#include <step_B.h>


void create_AQ_pairs(const Params P, mpz_t* Ans, mpz_t* Qns, unsigned f,
					 mpz_t* exp_vects, mpz_t* hist_vects,
					 const mpz_t* factor_base, const unsigned B) {
/*
	TODO: docstring me

	:param P: Set of parameters for the problem (see step_A.h).
	:param f: Numb. of factored Qns.
	:param B: Cardinality of factor_base.
*/

	/****************************
	* Declarations, allocations *
	****************************/
	
	unsigned* Qn_odd_pows;  // odd_power
	unsigned n_Qn_odd_pows;  // s_odd_power
	//		/\
	// @Margot: j'ai mis n et non s car à mon avis s doit indiquer la
	// taille du tableau au sens de C et ici le tableau Qn_odd_pows
	// n'est pas de taille n_Qn_odd_pows à proprement parler, il est de
	// taille B.
	unsigned* reduced_fb_indexes;  // prime_table
	unsigned s_reduced_fb_indexes;  // s_prime_table

	Qn_odd_pows = (unsigned*)malloc(B * sizeof(unsigned));
	reduced_fb_indexes = (unsigned*)malloc(B * sizeof(unsigned));

	/*******
	* Free *
	*******/

	free(Qn_odd_pows);
	free(reduced_fb_indexes);
}
