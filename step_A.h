/*
	step_A.h
*/

#ifndef STEP_A
#define STEP_A


struct Params {
	mpz_t	 N;				  // Integer to factor
	unsigned k;				  // Integer k used in sqrt(kN)
	mpz_t    n_lim;			  // Size of the cont. frac. expension of sqrt(kN)
	unsigned n_want_factQns;  // Number of wanted factored Qns
};

typedef struct Params Params
// @Funk Voir https://c.programmingpedia.net/fr/tutorial/2681/typedef

#endif
