#ifndef BARRET_H
#define BARRET_H

#include <gmp.h>
#include <stdint.h>



void limb_sum_mod64p1_64m1(mpz_t N, uint64_t* sum64min1, uint64_t* sum64min2);

void barrett_reduction_OURS(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n);

void barrett_reduction_REF(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n);

void barret_multiplication(mpz_t result, mpz_t a, mpz_t b, mpz_t q, size_t k);

void barret_multiplication_radix(mpz_t result, const mpz_t X, const mpz_t Y, const mpz_t M, unsigned long m, unsigned long alpha, unsigned long beta);

#endif //BARRET_H
