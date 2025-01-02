#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define WORD_SIZE 32
#define SIZE_U 4096
#define SIZE_N 2048

void compute_sums(mpz_t N, int *sum1, int *sum2, int *sum3) {
    // Variables to store intermediate sums
    int sum_temp0 = 0, sum_temp1 = 0, sum_temp2 = 0, sum_temp3 = 0;

    // Convert N to a string for digit processing
    char *digits = mpz_get_str(NULL, 10, N);

    // Process digits from least significant to most significant
    size_t len = strlen(digits);
    for (size_t i = 0; i < len; i++) {
        int digit = digits[len - i - 1] - '0'; // Get the numeric value of the digit

        if (i % 4 == 0) {
            sum_temp0 += digit;
        } else if (i % 4 == 1) {
            sum_temp1 += digit;
        } else if (i % 4 == 2) {
            sum_temp2 += digit;
        } else if (i % 4 == 3) {
            sum_temp3 += digit;
        }
    }

    // Free the string allocated by mpz_get_str
    free(digits);

    // Compute final sums
    *sum1 = (sum_temp0 + sum_temp1 + sum_temp2 + sum_temp3) % 9;                  // Mod 9
    *sum2 = (sum_temp0 - sum_temp1 + sum_temp2 - sum_temp3) % 11;                 // Mod 11
    *sum3 = (sum_temp0 + (10 * sum_temp1) - sum_temp2 - (10 * sum_temp3)) % 101;  // Mod 101
    // Ensure sums are positive
    if (*sum1 < 0) *sum1 += 9;
    if (*sum2 < 0) *sum2 += 11;
    if (*sum3 < 0) *sum3 += 101;
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


int main() {
    mpz_t N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett;
    unsigned int n, length_N, counter_result = 0, iteration = 10000;
    int sum_N_1, sum_N_2, sum_N_3, sum_q_1, sum_q_2, sum_q_3, sum_u_1, sum_u_2, sum_u_3, sum_r_1, sum_r_2, sum_r_3;
    int Zero_1, Zero_2, Zero_3;
    int counter1 = 0, counter2 = 0, counter3 = 0, counter_loop = 0, counter_combined = 0;
    int fault_happened_loop = 0;

    // Initialize GMP integers
    mpz_inits(N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, NULL);

    // Initialize random state for GMP
    gmp_randstate_t state;
    gmp_randinit_mt(state); // Mersenne Twister random generator
    gmp_randseed_ui(state, time(NULL)); // Seed with current time


    // Iterate for simulations
    for (unsigned int i = 0; i < iteration; i++) {
        // Correct values

        // Generate a random N
        mpz_urandomb(N, state, SIZE_N);

        // Generate a random u
        mpz_urandomb(u, state, SIZE_U);

        // Calculate length_N and n
        length_N = mpz_sizeinbase(N, 2); // Number of bits in N
        n = (length_N + WORD_SIZE - 1) / WORD_SIZE;

        // Precompute R, b^(n-1), b^(n+1)
        mpz_ui_pow_ui(b_n_minus_1, 2, (n - 1) * WORD_SIZE);
        mpz_ui_pow_ui(b_n_plus_1, 2, (n + 1) * WORD_SIZE);
        mpz_ui_pow_ui(R, 2, 2 * n * WORD_SIZE);
        mpz_fdiv_q(R, R, N); // R = floor(b^(2n) / N)

/**************************************************************************************************************************/
        // Perform Barrett reduction
        barrett_reduction_OURS(r_barrett, q_barrett, &fault_happened_loop, u, N, R, b_n_minus_1, b_n_plus_1, n);
        counter_loop = counter_loop + fault_happened_loop;
        // printf("fault_happened_loop: %d\n", fault_happened_loop);


        compute_sums(N, &sum_N_1, &sum_N_2, &sum_N_3);
        compute_sums(u, &sum_u_1, &sum_u_2, &sum_u_3);
        compute_sums(r_barrett, &sum_r_1, &sum_r_2, &sum_r_3);
        compute_sums(q_barrett, &sum_q_1, &sum_q_2, &sum_q_3);


        Zero_1 = (sum_u_1) - (sum_N_1 * sum_q_1) - sum_r_1;

        Zero_2 = (sum_u_2) - (sum_N_2 * sum_q_2) - sum_r_2;

        Zero_3 = (sum_u_3) - (sum_N_3 * sum_q_3) - sum_r_3;


        if (Zero_1 % 9 == 0 && Zero_2 % 11 == 0 && Zero_3 % 101 == 0 && fault_happened_loop == 0) {
            counter_combined++;
        }
    }

/**************************************************************************************************************************/

        // barrett_reduction_REF(r_barrett, u, N, R, b_n_minus_1, b_n_plus_1, n);


    // Print results
    printf("No fault happened is %u.\n", counter_result); // only for checking it works- must delete for overhead
    // printf("Not detected from sum1 is %u.\n", counter1);
    // printf("Not detected from sum2 is %u.\n", counter2);
    // printf("Not detected from sum3 is %u.\n", counter3);
    printf("fault happened but not detected from combined is %u.\n", counter_combined);

    // Clear GMP integers
    mpz_clears(N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, NULL);

    return 0;
}
