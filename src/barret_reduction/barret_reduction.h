#ifndef BARRET_REDUCTION_H
#define BARRET_REDUCTION_H

#include <gmp.h>

/**
 * @brief Barret Reduction with countermeasure added
 * @param r Parameter r
 * @param q Parameter q
 * @param fault_happened If fault happened, set to 1
 * @param u Parameter u
 * @param N Parameter N
 * @param R Parameter R
 * @param b_n_plus_1 Value
 * @param n Parameter n
 */
void BR_Robust(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N,
               const mpz_t R, const mpz_t b_n_plus_1, unsigned int n);

/**
 * @brief Original implementation of the barret reduction
 * @param r Parameter r
 * @param u Parameter u
 * @param N Parameter N
 * @param R Parameter R
 * @param b_n_plus_1 Value
 * @param n Parameter n
 */
void BR_Origin(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n);

#endif
