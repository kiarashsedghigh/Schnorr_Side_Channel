#include <gmp.h>
#include "size.h"

void BR_Robust(mpz_t r, mpz_t q, int *fault_happened, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1,
             unsigned int n) {
    mpz_t u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2;
    unsigned int subtraction_count = 0;

    // Initialize temporary variables
    mpz_inits(u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2, NULL);

    // Step 1: Compute q_hat = floor(floor(u / b^(n-1)) * R / b^(n+1))
    mpz_tdiv_q_2exp(u_high, u, (n - 1) * WORD_SIZE); // u_high = u >> (n-1)*word_size
    mpz_mul(q_hat, u_high, R); // q_hat = u_high * R
    mpz_tdiv_q_2exp(q_hat, q_hat, (n + 1) * WORD_SIZE); // q_hat = q_hat >> (n+1)*word_size

    // Step 2: Compute r = u mod b^(n+1) - (q_hat * N) mod b^(n+1)
    mpz_sub_ui(b_n_plus_2, b_n_plus_1, 1); // b_n_plus_1 - 1

    // r1 = u & (b_n_plus_1 - 1) (equivalent to u % b^(n+1))
    mpz_and(r1, u, b_n_plus_2);

    // r2 = (q_hat * N) & (b_n_plus_1 - 1) (equivalent to (q_hat * N) % b^(n+1))
    mpz_mul(N_mult, q_hat, N); // N_mult = q_hat * N
    mpz_and(r2, N_mult, b_n_plus_2); // r2 = (q_hat * N) & (b_n_plus_1 - 1)

    // r = r1 - r2
    mpz_sub(r, r1, r2);

    // Step 3: If r < 0, add b^(n+1)
    if (mpz_sgn(r) < 0)
        mpz_add(r, r, b_n_plus_2);


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

void BR_Origin(mpz_t r, const mpz_t u, const mpz_t N, const mpz_t R, const mpz_t b_n_plus_1, unsigned int n) {
    mpz_t u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2;

    // Initialize temporary variables
    mpz_inits(u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2, NULL);

    // Step 1: Compute q_hat = floor(floor(u / b^(n-1)) * R / b^(n+1))
    mpz_tdiv_q_2exp(u_high, u, (n - 1) * WORD_SIZE); // u_high = u >> (n-1)*word_size
    mpz_mul(q_hat, u_high, R); // q_hat = u_high * R
    mpz_tdiv_q_2exp(q_hat, q_hat, (n + 1) * WORD_SIZE); // q_hat = q_hat >> (n+1)*word_size

    // Step 2: Compute r = u mod b^(n+1) - (q_hat * N) mod b^(n+1)
    // Instead of using mpz_mod, use bitwise AND
    mpz_sub_ui(b_n_plus_2, b_n_plus_1, 1); // b_n_plus_1 - 1

    // r1 = u & (b_n_plus_1 - 1) (equivalent to u % b^(n+1))
    mpz_and(r1, u, b_n_plus_2);

    // r2 = (q_hat * N) & (b_n_plus_1 - 1) (equivalent to (q_hat * N) % b^(n+1))
    mpz_mul(N_mult, q_hat, N); // N_mult = q_hat * N
    mpz_and(r2, N_mult, b_n_plus_2); // r2 = (q_hat * N) & (b_n_plus_1 - 1)

    // r = r1 - r2
    mpz_sub(r, r1, r2);

    // Step 3: If r < 0, add b^(n+1)
    if (mpz_sgn(r) < 0)
        mpz_add(r, r, b_n_plus_2);


    // Step 4: While r >= N, subtract N
    while (mpz_cmp(r, N) >= 0)
        mpz_sub(r, r, N);


    mpz_clears(u_high, q_hat, r1, r2, temp, N_mult, NULL);
}


