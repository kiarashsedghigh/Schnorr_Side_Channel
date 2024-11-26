#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <papi.h>
#include "src/barret.h"
#include "src/size.h"

int main() {
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


    // Initialize the PAPI library
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI library initialization error!\n");
        exit(1);
    }
    int EventSet = PAPI_NULL;
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr, "Error creating PAPI event set\n");
        exit(1);
    }

    // Add events for total cycles and total instructions
    if (PAPI_add_event(EventSet, PAPI_TOT_CYC) != PAPI_OK) {
        fprintf(stderr, "Error adding PAPI_TOT_CYC\n");
        exit(1);
    }
    if (PAPI_add_event(EventSet, PAPI_TOT_INS) != PAPI_OK) {
        fprintf(stderr, "Error adding PAPI_TOT_INS\n");
        exit(1);
    }

    /* Variables holding the cycle/instruction counts for single/total interation(s) */
    long long benchmark_results_total_iterations[2] = {0,0};
    long long benchmark_results_single_iteration[2] = {0,0};

    // Iterate for simulations
    for (unsigned int i = 0; i < iteration; i++) {

        // Generate a random N
        mpz_urandomb(N, state, SIZE_N);

        // Generate a random u
        mpz_urandomb(u, state, SIZE_U);



/**************************************************************************************************************************/
        // Start the event set
        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "Error starting PAPI\n");
            exit(1);
        }

        // Calculate length_N and n
        unsigned int length_N = mpz_sizeinbase(N, 2); // Number of bits in N
        n = (length_N + WORD_SIZE - 1) / WORD_SIZE;

        // Precompute R, b^(n-1), b^(n+1)
        mpz_ui_pow_ui(b_n_minus_1, 2, (n - 1) * WORD_SIZE);
        mpz_ui_pow_ui(b_n_plus_1, 2, (n + 1) * WORD_SIZE);
        mpz_ui_pow_ui(R, 2, 2 * n * WORD_SIZE);
        mpz_fdiv_q(R, R, N); // R = floor(b^(2n) / N)

        // Perform Barrett reduction

        barrett_reduction_OURS(r_barrett, q_barrett, &fault_happened_loop, u, N, R, b_n_minus_1, b_n_plus_1, n);
        // barrett_reduction_REF(r_barrett, u, N, R, b_n_minus_1, b_n_plus_1, n);

        // Adding is done in the function
        counter_loop = counter_loop + fault_happened_loop;


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
