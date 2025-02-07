#ifndef BARRET_MULTIPLICATION_H
#define BARRET_MULTIPLICATION_H

#include <gmp.h>

/**
 * @brief Encode the given numbers based on the method presented in the paper
 * @param X First input number
 * @param Y Second input number
 * @param X_encoded Encoding of the first number
 * @param Y_encoded Encoding of the second number
 * @param M Modulo
 * @param shift Number of bits to shift the second number with
*/
void BM_encode_inputs(const mpz_t X, const mpz_t Y, mpz_t X_encoded, mpz_t Y_encoded, mpz_t M, int shift);


/**
 * @brief Barret Multiplication main function
 * @param result Result of (X * Y) % M
 * @param partial_computation_result Result of intermediate steps of the computation (Z_i value)
 * @param X First input number
 * @param Y Second input number
 * @param M Modulo
 * @param m Parameter 'm' of the Barret Multiplication
 * @param alpha Parameter 'a' of the Barret Multiplication
 * @param beta Parameter 'b' of the Barret Multiplication
 * @param recompute_count Number of steps to be recomputed
 */
void BM_radix(mpz_t result, mpz_t partial_computation_result, const mpz_t X, const mpz_t Y, const mpz_t M,
              unsigned long m,
              unsigned long alpha, long beta, int recompute_count);


/**
 * @brief Barret Multiplication recomputation function
 * @param partial_computation_result Result of intermediate steps of the computation (Z_i value)
 * @param X First input number
 * @param Y Second input number
 * @param M Modulo
 * @param m Parameter 'm' of the Barret Multiplication
 * @param alpha Parameter 'a' of the Barret Multiplication
 * @param beta Parameter 'b' of the Barret Multiplication
 * @param recompute_count Number of steps to be recomputed
 * @param Y_bit_length Length of the second input in bits
 */
void BM_radix_recomputation(mpz_t partial_computation_result, const mpz_t X, const mpz_t Y, const mpz_t M,
                            unsigned long m, unsigned long alpha, long beta, int recompute_count, int Y_bit_length);


/**
 * @brief Setup function to precompute input-independent values
 * @param M Modulo
 * @param m Parameter 'm' of the Barret Multiplication
 * @param alpha Parameter 'a' of the Barret Multiplication
 */
void setup_BM_radix(const mpz_t M, unsigned long m, unsigned long alpha);

/**
 * @brief Deallocated memories allocated during the setup phase
 */
void setup_clear_BM_radix();


#endif
