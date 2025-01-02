#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <papi.h>
#include "src/barret.h"
#include "src/size.h"
#include <stdint.h>



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

    mpz_t N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett;
    unsigned int n, counter_result = 0, iteration = 10000;
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


    __uint128_t mod_64plus1 = ((__uint128_t)1 << 64) + 1; // Use 128-bit integer
    __uint128_t mod_64minus1 = 18446744073709551615; // Use 128-bit integer
    __uint128_t mod_64minus2 = 18446744073709551614; // Use 128-bit integer

    init_r();
    /* Variables holding the cycle/instruction counts for single/total interation(s) */
    long long benchmark_results_total_iterations[2] = {0,0};
    long long benchmark_results_single_iteration[2] = {0,0};


    // Iterate for simulations
    for (unsigned int i = 0; i < iteration; i++) {

        // Generate a random N
        mpz_urandomb(N, state, SIZE_N);

        // Generate a random u
        mpz_urandomb(u, state, SIZE_U);
        mpz_mod(r_correct, u, N);  // r_correct = u % N
        mpz_fdiv_q(q_correct, u, N); // q_correct = u // N


        uint64_t sum1_N = sum_eff_mod64min1(N);
        uint64_t sum2_N = sum_eff_mod64min2(N);



/**************************************************************************************************************************/

        // Calculate length_N and n
        unsigned int length_N = mpz_sizeinbase(N, 2); // Number of bits in N
        n = (length_N + WORD_SIZE - 1) / WORD_SIZE;

        // Precompute R, b^(n-1), b^(n+1)
        // mpz_ui_pow_ui(b_n_minus_1, 2, (n - 1) * WORD_SIZE);
        mpz_ui_pow_ui(b_n_plus_1, 2, (n + 1) * WORD_SIZE);
        mpz_ui_pow_ui(R, 2, 2 * n * WORD_SIZE);
        mpz_fdiv_q(R, R, N); // R = floor(b^(2n) / N)

        mpz_t diffur;
        mpz_init(diffur);

        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "Error starting PAPI\n");
            exit(1);
        }


        barrett_reduction_OURS(r_barrett, q_barrett, &fault_happened_loop, u, N, R, b_n_plus_1, n);


        //After Barret

        mpz_sub(diffur, u, r_barrett);

        uint64_t sum1_ur = sum_eff_mod64min1(diffur);
        uint64_t sum1_qbarret = sum_eff_mod64min1(q_barrett);


        uint64_t sum2_ur = sum_eff_mod64min2(diffur);
        uint64_t sum2_qbarret = sum_eff_mod64min2(q_barrett);

        __uint128_t mul1_sqsn = (__uint128_t)sum1_qbarret * (__uint128_t)sum1_N;
        __uint128_t zero1 = ((__uint128_t)(mul1_sqsn) - (__uint128_t)sum1_ur) % (mod_64minus1);


        __uint128_t mul2_sqsn = (__uint128_t)sum2_qbarret * (__uint128_t)sum2_N;
        __uint128_t zero2 =  ((__uint128_t)(mul2_sqsn) - (__uint128_t)sum2_ur) % (mod_64minus2);


        if (zero1 == 0 && zero2 == 0 && fault_happened_loop == 0) {
            counter_combined++;
        }

        if (mpz_cmp(r_barrett, r_correct) == 0) {
            counter_result++;
        }

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

    // Display results
    printf("Total cycles: %lld\n", benchmark_results_total_iterations[0] / iteration);
    printf("Total instructions: %lld\n", benchmark_results_total_iterations[1] / iteration);

    /**************************************************************************************************************************/
    printf("No fault happened is %u.\n", counter_result); // only for checking it works- must delete for overhead
    printf("fault happened but not detected from combined is %u.\n", counter_combined);

    // Clear GMP integers
    mpz_clears(N, u, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, NULL);
    return 0;


}