#ifndef BARRET_H
#define BARRET_H

#include <gmp.h>
#include <stdint.h>



void sum_eff_mod64twomodes(mpz_t N, uint64_t* sum64min1, uint64_t* sum64min2);

void barrett_reduction_OURS(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n);

void barrett_reduction_REF(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n);



#endif //BARRET_H
