#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "size.h"
#include <stdint.h>

__uint128_t mod_64minus1 = ((__uint128_t)1 << 64) - 1; // Use 128-bit integer
__uint128_t mod_64minus2 = ((__uint128_t)1 << 64) - 2; // Use 128-bit integer

void limb_sum_mod64p1_64m1(mpz_t N, uint64_t* sum64min1, uint64_t* sum64min2) {
    const mp_size_t limb_cnt = mpz_size(N);  // Number of limbs in the mpz_t number
    const mp_limb_t *limbs = mpz_limbs_read(N);  // Get the limbs of the number
    __uint128_t temp_sum_64min1 = 0;
    __uint128_t temp_sum_64min2 = 0;

    for (int i = 0; i < limb_cnt; i++) {
        // 2^64 - 1
        temp_sum_64min1 += (__uint128_t)limbs[i];

        // 2^64 - 2
        temp_sum_64min2 += (((__uint128_t)((__uint128_t)limbs[i] << i)));
    }

    *sum64min1 = temp_sum_64min1 % mod_64minus1;
    *sum64min2 = temp_sum_64min2 % mod_64minus2;
}


void barrett_reduction_OURS(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n) {

    mpz_t u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2;
    unsigned int subtraction_count = 0;

    // Initialize temporary variables
    mpz_inits(u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2, NULL);

    // Step 1: Compute q_hat = floor(floor(u / b^(n-1)) * R / b^(n+1))
    mpz_tdiv_q_2exp(u_high, u, (n - 1) * WORD_SIZE); // u_high = u >> (n-1)*word_size
    mpz_mul(q_hat, u_high, R);                       // q_hat = u_high * R
    mpz_tdiv_q_2exp(q_hat, q_hat, (n + 1) * WORD_SIZE); // q_hat = q_hat >> (n+1)*word_size

    // Step 2: Compute r = u mod b^(n+1) - (q_hat * N) mod b^(n+1)
    // Instead of using mpz_mod, use bitwise AND
    mpz_sub_ui(b_n_plus_2, b_n_plus_1, 1); // b_n_plus_1 - 1

    // r1 = u & (b_n_plus_1 - 1) (equivalent to u % b^(n+1))
    mpz_and(r1, u, b_n_plus_2);

    // r2 = (q_hat * N) & (b_n_plus_1 - 1) (equivalent to (q_hat * N) % b^(n+1))
    mpz_mul(N_mult, q_hat, N);  // N_mult = q_hat * N
    mpz_and(r2, N_mult, b_n_plus_2);  // r2 = (q_hat * N) & (b_n_plus_1 - 1)

    // r = r1 - r2
    mpz_sub(r, r1, r2);


    // Step 3: If r < 0, add b^(n+1)
    if (mpz_sgn(r) < 0) {
        mpz_add(r, r, b_n_plus_2);
    }

    // Step 4: While r >= N, subtract N
    while (mpz_cmp(r, N) >= 0) {
        mpz_sub(r, r, N);
        subtraction_count++;
        if (subtraction_count > 2) {
            *fault_happened = 1;
            break;
        }
    }

    // Correct q = q_hat + subtraction_count
    mpz_add_ui(q, q_hat, subtraction_count);

    // Clear temporary variables
    mpz_clears(u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2, NULL);
}
//
void barrett_reduction_REF(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n) {
    mpz_t u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2;

    // Initialize temporary variables
    mpz_inits(u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2, NULL);

    // Step 1: Compute q_hat = floor(floor(u / b^(n-1)) * R / b^(n+1))
    mpz_tdiv_q_2exp(u_high, u, (n - 1) * WORD_SIZE); // u_high = u >> (n-1)*word_size
    mpz_mul(q_hat, u_high, R);                       // q_hat = u_high * R
    mpz_tdiv_q_2exp(q_hat, q_hat, (n + 1) * WORD_SIZE); // q_hat = q_hat >> (n+1)*word_size

    // Step 2: Compute r = u mod b^(n+1) - (q_hat * N) mod b^(n+1)
    // Instead of using mpz_mod, use bitwise AND
    mpz_sub_ui(b_n_plus_2, b_n_plus_1, 1); // b_n_plus_1 - 1

    // r1 = u & (b_n_plus_1 - 1) (equivalent to u % b^(n+1))
    mpz_and(r1, u, b_n_plus_2);

    // r2 = (q_hat * N) & (b_n_plus_1 - 1) (equivalent to (q_hat * N) % b^(n+1))
    mpz_mul(N_mult, q_hat, N);  // N_mult = q_hat * N
    mpz_and(r2, N_mult, b_n_plus_2);  // r2 = (q_hat * N) & (b_n_plus_1 - 1)

    // r = r1 - r2
    mpz_sub(r, r1, r2);

    // Step 3: If r < 0, add b^(n+1)
    if (mpz_sgn(r) < 0) {
        mpz_add(r, r, b_n_plus_2);
    }

    // Step 4: While r >= N, subtract N
    while (mpz_cmp(r, N) >= 0) {
        mpz_sub(r, r, N);
    }

    mpz_clears(u_high, q_hat, r1, r2, temp, N_mult, NULL);
}



void barret_multiplication_radix(mpz_t result, const mpz_t X, const mpz_t Y, const mpz_t M, unsigned long m, unsigned long alpha, unsigned long beta) {
    mpz_t r, mu, T, q, temp;
    mpz_inits(r, mu, T, q, temp, NULL);

    // Step 1: Precompute r = 2^m and mu = floor(2^(N+alpha) / M)
    mpz_ui_pow_ui(r, 2, m);
    size_t N = mpz_sizeinbase(M, 2); // N = bit length of M
    mpz_ui_pow_ui(temp, 2, N + alpha);
    mpz_fdiv_q(mu, temp, M);

    // Step 2: Initialize Z = 0
    mpz_t Z;
    mpz_init_set_ui(Z, 0);

    // Decompose Y into nm words
    size_t nm = (N + m - 1) / m; // ceiling(N / m)
    mpz_t *Y_words = malloc(nm * sizeof(mpz_t));
    for (size_t i = 0; i < nm; i++) {
        mpz_init(Y_words[i]);
        mpz_fdiv_r_2exp(Y_words[i], Y, m * (i + 1));
        mpz_fdiv_q_2exp(Y_words[i], Y_words[i], m * i);
    }

    // Step 3: Main loop
    for (ssize_t i = nm - 1; i >= 0; i--) {
        mpz_mul(T, Z, r);          // T = Z(i+1) * r
        mpz_addmul(T, X, Y_words[i]); // T += X * Y_i
        mpz_mul(temp, T, mu);      // temp = T * mu
        mpz_fdiv_q_2exp(q, temp, N + beta); // q = floor(temp / 2^(N+beta))
        mpz_mul(temp, q, M);       // temp = q * M
        mpz_sub(Z, T, temp);       // Z = T - temp
    }

    // Step 4: Final reduction
    if (mpz_cmp(Z, M) >= 0) {
        mpz_sub(Z, Z, M);
    }

    // Set result
    mpz_set(result, Z);

    // Cleanup
    mpz_clears(r, mu, T, q, temp, Z, NULL);
    for (size_t i = 0; i < nm; i++) {
        mpz_clear(Y_words[i]);
    }
    free(Y_words);
}



void barret_multiplication(mpz_t result, mpz_t a, mpz_t b, mpz_t q, size_t k) {

    // Initialize variables
    mpz_t t1, t2, t3, t4, m, two_k;
    mpz_inits(t1, t2, t3, t4, m, two_k, NULL);

    // Compute k = ceil(log2(q))
    // k = mpz_sizeinbase(q, 2);

    // Compute m = floor(2^(2k) / q)
    mpz_ui_pow_ui(two_k, 2, 2 * k); // two_k = 2^(2k)
    mpz_fdiv_q(m, two_k, q);        // m = floor(2^(2k) / q)

    // Compute t1 = a * b
    mpz_mul(t1, a, b);

    // Compute t2 = floor(t1 * m / 2^(2k))
    mpz_mul(t2, t1, m);         // t2 = t1 * m
    mpz_fdiv_q_2exp(t2, t2, 2 * k); // t2 = floor(t2 / 2^(2k))

    // Compute t3 = t2 * q
    mpz_mul(t3, t2, q);

    // Compute t4 = t1 - t3
    mpz_sub(t4, t1, t3);

    // If t4 >= q, compute z = t4 - q
    if (mpz_cmp(t4, q) >= 0) {
        mpz_sub(result, t4, q);
    } else {
        mpz_set(result, t4);
    }

    // Clear variables
    mpz_clears(t1, t2, t3, t4, m, two_k, NULL);
}


// Precompute mu = floor(2^k / m) for a given modulus m
void compute_shoup_mu(mpz_t mu, const mpz_t m, unsigned int k) {
    mpz_t base;
    mpz_init(base);
    mpz_ui_pow_ui(base, 2, k);  // base = 2^k
    mpz_fdiv_q(mu, base, m);    // mu = floor(2^k / m)
    mpz_clear(base);
}

// Shoup modular multiplication
void shoup_multiply(mpz_t result, const mpz_t a, const mpz_t b, const mpz_t m, const mpz_t mu, unsigned int k) {
    mpz_t r, q, temp;
    mpz_inits(r, q, temp, NULL);

    // Step 1: r = a * b
    mpz_mul(r, a, b);

    // Step 2: q = floor(r * mu / 2^k)
    mpz_mul(temp, r, mu);        // temp = r * mu
    mpz_tdiv_q_2exp(q, temp, k); // q = temp / 2^k (right shift by k bits)

    // Step 3: result = r - q * m
    mpz_mul(temp, q, m);         // temp = q * m
    mpz_sub(result, r, temp);    // result = r - temp

    // Step 4: Normalize result
    if (mpz_cmp(result, m) >= 0) {
        mpz_sub(result, result, m); // result -= m if result >= m
    }
    if (mpz_sgn(result) < 0) {
        mpz_add(result, result, m); // result += m if result < 0
    }

    mpz_clears(r, q, temp, NULL);
}

