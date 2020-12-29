/*
	step_B.c
*/

#include <gmp.h>

#include <step_A.h>
#include <step_B.h>

#define LIM (


void create_AQ_pairs(const Params P, mpz_t* Ans, mpz_t* Qns, unsigned* f,
					 mpz_t* exp_vects, mpz_t* hist_vects,
					 const mpz_t* factor_base, const unsigned B) {
/*
	TODO: Document me

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
	unsigned n;

	Qn_odd_pows = (unsigned*)malloc(B * sizeof(unsigned));
	reduced_fb_indexes = (unsigned*)malloc(B * sizeof(unsigned));
	n = 0;

	// TODO: Expand sqrt(kN), step 2) chez Margot
	// TODO: Limite dynamique, step 3) chez Margot
	
	while (n < P.n_lim
			&& &f < P.n_want_fact_Qns != 0 ? P.n_want_fact_Qns : dyn_lim
			&& q_last != 1) {
		// TODO: Définir dyn_lim (changer de nom ?) et q_last

		// TODO: Calculer An-1 et Qn, step a) chez Margot
		// On les note Anm1 et Qn dans la suite du programme

		Ans[f] = Anm1;
		Qns[f] = Qn;

		// TODO: Appeler les add_*

		&f ++;
		n ++;
	}
			

	/*******
	* Free *
	*******/

	free(Qn_odd_pows);
	free(reduced_fb_indexes);
}
