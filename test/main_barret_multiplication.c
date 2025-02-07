#include "barret_lib.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define ITERATIONS 100000

int main() {
    int EventSet = PAPI_NULL;

    /* Initialize the PAPI library */
    if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
        fprintf(stderr, "PAPI initialization error.\n");
        return -1;
    }

    /* Create an Event Set */
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        fprintf(stderr, "Error creating event set.\n");
        return -1;
    }

    /* Add PAPI_TOT_CYC to the Event Set */
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

    /* Initialize GMP integers */
    mpz_t N, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, ur_diff;
    mpz_inits(N, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, ur_diff, NULL);

    mpz_t X, X_encoded, Y, Y_encoded, M, result, result_partial, result_partial2;
    mpz_inits(X, Y, X_encoded, Y_encoded, M, result, result_partial, result_partial2, NULL);


    /* Initialize random state for GMP */
    gmp_randstate_t state;
    gmp_randinit_mt(state); // Mersenne Twister random generator
    gmp_randseed_ui(state, time(NULL)); // Seed with current time

    unsigned int counter_result = 0;


    for (unsigned int i = 0; i < ITERATIONS; i++) {
        /* Generate a random u */
        mpz_urandomb(X, state, SIZE_U);
        mpz_urandomb(Y, state, SIZE_U);
        mpz_urandomb(M, state, SIZE_U);

        /* Perform Barrett Modular Multiplication */
        const unsigned long m = 64; // word size
        long alpha = m + 40, beta = -20; // parameters
        const int recompute_iter = 5;
        const int shift = 2;

        setup_BM_radix(M, m, alpha);

        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "Error starting PAPI\n");
            exit(1);
        }

        BM_radix(result, result_partial, X, Y, M, m, alpha, beta, recompute_iter);

        BM_encode_inputs(X, Y, X_encoded, Y_encoded, M, shift);

        BM_radix_recomputation(result_partial2, X_encoded, Y_encoded, M, m, alpha, beta, recompute_iter,
                               mpz_sizeinbase(Y, 2) + shift);

        if (mpz_cmp(result_partial, result_partial2) == 0)
            counter_result++;
        else {
            printf("Recomputation Fault Detected \n");
            exit(0);
        }

        if (PAPI_stop(EventSet, benchmark_results_single_iteration) != PAPI_OK) {
            fprintf(stderr, "Error stopping PAPI\n");
            exit(1);
        }

        benchmark_results_total_iterations[0] += benchmark_results_single_iteration[0];
        benchmark_results_total_iterations[1] += benchmark_results_single_iteration[1];

        setup_clear_BM_radix();
    }
    /* Cleanup */
    mpz_clears(X, Y, M, result, result_partial, result_partial2, NULL);


    /* Cleanup */
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
    mpz_clears(N, R, b_n_minus_1, b_n_plus_1, r_correct, q_correct, r_barrett, q_barrett, NULL);
    return 0;
}
