#ifndef BARRET_H
#define BARRET_H

#include <gmp.h>



void compute_sums(mpz_t N, int *sum1, int *sum2, int *sum3);

void barrett_reduction_OURS(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_minus_1, const mpz_t b_n_plus_1, unsigned int n);

void barrett_reduction_REF(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_minus_1, const mpz_t b_n_plus_1, unsigned int n);


#endif //BARRET_H
