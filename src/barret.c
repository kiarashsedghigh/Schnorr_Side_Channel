#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "size.h"
#include <stdint.h>


uint64_t barrett_reduction_64bit_mod_10(uint64_t x) {
    const uint64_t mu = (1ULL << 64) / 10;  // mu = floor(2^64 / 10)

    const uint64_t q = (x * mu) >> 64;      // q = floor(x * mu / 2^64)
    uint64_t r = x - q * 10;          // r = x - q * 10

    if (r >= 10)
        r -= 10;

    return r;
}

void compute_sums(mpz_t N, int *sum1, int *sum2, int *sum3) {

    int sum_temps[4] = {0, 0, 0, 0};  // Array to store intermediate sums

    mp_size_t limb_cnt = mpz_size(N);  // Number of limbs in the mpz_t number

    mp_limb_t *limbs = mpz_limbs_read(N);  // Get the limbs of the number

    size_t digit_count = 0;

    // Process each limb in the number (from least significant to most significant)
    for (mp_size_t i = 0; i < limb_cnt; i++) {
        mp_limb_t limb = limbs[i];  // Get the current limb

        // Extract digits from the limb
        while (limb > 0) {
            const int digit = (int) barrett_reduction_64bit_mod_10(limb);

            sum_temps[digit_count % 4] += digit;  // Add the digit to the corresponding sum bucket

            // if (digit_count == 4)
            //     digit_count = 0;
            limb /= 10;  // Divide the limb by 10 to process the next digit

            digit_count++;
        }
    }
    *sum1 = (sum_temps[0] + sum_temps[1] + sum_temps[2] + sum_temps[3] + 9) % 9;                  // Mod 9
    *sum2 = (sum_temps[0] - sum_temps[1] + sum_temps[2] - sum_temps[3] + 11) % 11;                 // Mod 11
    *sum3 = (sum_temps[0] + (10 * sum_temps[1]) - sum_temps[2] - (10 * sum_temps[3]) + 101) % 101;  // Mod 101
}


void barrett_reduction_OURS(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_minus_1, const mpz_t b_n_plus_1, unsigned int n) {

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

void barrett_reduction_REF(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_minus_1, const mpz_t b_n_plus_1, unsigned int n) {
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
