#include "barret_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ITERATIONS 100000


int main() {
    int EventSet = PAPI_NULL;

    // Initialize the PAPI library
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI initialization error.\n");
        return -1;
    }

    // Create an Event Set
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr, "Error creating event set.\n");
        return -1;
    }

    // Add PAPI_TOT_CYC to the Event Set
    if (PAPI_add_event(EventSet, PAPI_TOT_CYC) != PAPI_OK) {
        fprintf(stderr, "Error adding PAPI_TOT_CYC.\n");
        return -1;
    }

    if (PAPI_add_event(EventSet, PAPI_TOT_INS) != PAPI_OK) {
        fprintf(stderr, "Error adding PAPI_TOT_INS\n");
        exit(1);
    }

    /* Variables holding the cycle/instruction counts for single/total interation(s) */
    long long benchmark_results_single_iteration[2] = {0, 0};
    long long benchmark_results_total_iterations[2] = {0, 0};


    // Initialize GMP integers
    mpz_t N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, ur_diff;
    mpz_inits(N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, ur_diff, NULL);

    // Initialize random state for GMP
    gmp_randstate_t state;
    gmp_randinit_mt(state); // Mersenne Twister random generator
    gmp_randseed_ui(state, time(NULL)); // Seed with current time

    unsigned int counter_result = 0;
    int counter_combined = 0;
    int fault_happened_loop = 0;


    for (unsigned int i = 0; i < ITERATIONS; i++) {
        /* Generate random N and u and compute q,r as N = u.q + r */
        // Generate a random N
        mpz_urandomb(N, state, SIZE_N);

        // Calculate number of words of number N
        const unsigned int number_word_count = (mpz_sizeinbase(N, 2) + WORD_SIZE - 1) / WORD_SIZE;

        // Generate a random u
        mpz_urandomb(u, state, SIZE_U);

        // Compute the real/correct value of q and r
        mpz_mod(r_correct, u, N);
        mpz_fdiv_q(q_correct, u, N);

        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "Error starting PAPI\n");
            exit(1);
        }

        // Precompute R, b^(n-1), b^(n+1)
        mpz_ui_pow_ui(b_n_plus_1, 2, (number_word_count + 1) * WORD_SIZE);
        mpz_ui_pow_ui(R, 2, 2 * number_word_count * WORD_SIZE);
        mpz_fdiv_q(R, R, N);

        /* Precompute barret ..... for already known and fixed N */
        uint64_t sum1_N, sum2_N;
        limb_sum_mod64m1_64m2(N, &sum1_N, &sum2_N);


        /* Calling Original Barret Reduction */
        // BR_Origin(r_barrett, u, N, R, b_n_plus_1, number_word_count);

        // /* Call our barret reduction implementation */
        // BR_Robust(r_barrett, q_barrett, &fault_happened_loop, u, N, R, b_n_plus_1, number_word_count);
        //
        // mpz_sub(ur_diff, u, r_barrett);
        //
        // uint64_t sum1_ur, sum2_ur;
        // uint64_t sum1_qbarret, sum2_qbarret;
        //
        //
        // limb_sum_mod64m1_64m2(ur_diff, &sum1_ur, &sum2_ur);
        // limb_sum_mod64m1_64m2(q_barrett, &sum1_qbarret, &sum2_qbarret);
        //
        // const __uint128_t mul1_sqsn = (__uint128_t) sum1_qbarret * (__uint128_t) sum1_N;
        // const __uint128_t zero1 = ((mul1_sqsn) - (__uint128_t) sum1_ur) % mod_64minus1;
        //
        // const __uint128_t mul2_sqsn = (__uint128_t) sum2_qbarret * (__uint128_t) sum2_N;
        // const __uint128_t zero2 = ((mul2_sqsn) - (__uint128_t) sum2_ur) % mod_64minus2;
        //
        // if (zero1 == 0 && zero2 == 0 && fault_happened_loop == 0)
        //     counter_combined++;

        if (PAPI_stop(EventSet, benchmark_results_single_iteration) != PAPI_OK) {
            fprintf(stderr, "Error stopping PAPI\n");
            exit(1);
        }

        benchmark_results_total_iterations[0] += benchmark_results_single_iteration[0];
        benchmark_results_total_iterations[1] += benchmark_results_single_iteration[1];
    }

    // Clean up
    if (PAPI_cleanup_eventset(EventSet) != PAPI_OK) {
        fprintf(stderr, "Error cleaning up event set\n");
        exit(1);
    }
    if (PAPI_destroy_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr, "Error destroying event set\n");
        exit(1);
    }
    PAPI_shutdown();


    /**************************************************************************************************************************/
    printf("Total cycles: %lld\n", benchmark_results_total_iterations[0] / ITERATIONS);
    printf("Total instructions: %lld\n", benchmark_results_total_iterations[1] / ITERATIONS);
    printf("No fault happened is %u.\n", counter_result);
    printf("fault happened but not detected from combined is %u.\n", counter_combined);

    mpz_clears(N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, NULL);
    return 0;
}
