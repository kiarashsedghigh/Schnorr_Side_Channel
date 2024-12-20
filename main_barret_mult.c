#include <stdio.h>
#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <papi.h>
#include "src/barret.h"
#include "src/size.h"

#include <stdint.h>

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
    long long benchmark_results_single_iteration[2] = {0,0};
    long long benchmark_results_total_iterations[2] = {0,0};



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

        const __uint128_t mod_64minus2 = 18446744073709551614;
        const __uint128_t mod_64minus1 = 18446744073709551615;

        /* Generate random N and u and compute q,r as N = u.q + r */
        // Generate a random N
        mpz_urandomb(N, state, SIZE_N);

        // Calculate number of words of number N
        const unsigned int number_word_count = (mpz_sizeinbase(N, 2) + WORD_SIZE - 1) / WORD_SIZE;

        // Generate a random u
        mpz_urandomb(u, state, SIZE_U);


        if (PAPI_start(EventSet) != PAPI_OK) {
            fprintf(stderr, "Error starting PAPI\n");
            exit(1);
        }

        // mpz_t a, b, q, mu, result;
        // mpz_inits(a, b, q, mu, result, NULL);
        //
        // // Set values (example inputs)
        // mpz_set_str(a, "12345678901234567829837489273490", 10);  // a
        // mpz_set_str(b, "387654321018923758927359283765892639876543210", 10);  // b
        // mpz_set_str(q, "884467440735123509817209847109274109723019273987123987709551617", 10);  // m = 2^64 + 1
        //
        // barret_multiplication(result, a, b, q, mpz_sizeinbase(q, 2));
        //
        // // gmp_printf("%Zd\n", result);
        // exit(0);


        mpz_t X, Y, M, result;
        mpz_inits(X, Y, M, result, NULL);

        // Initialize values for X, Y, M
        mpz_set_str(X, "1231892323137891723981293861824618264123145678912345678912318923231378917239812938618246182641231456789123456789", 10);
        mpz_set_str(Y, "1281739718273891273111111111111111111111111111111118924981648623846263412817397182738912731111111111111111111111111111111189249816486238462634", 10);
        mpz_set_str(M, "928375892749237498723948723984729834798237582623984618947123719748923648926734928375892749237498723948723984729834798237582623984618947123719748923648926734", 10);


        // Perform Barrett Modular Multiplication
        unsigned long m = 16; // word size
        unsigned long alpha = 16, beta = 16; // parameters
        barret_multiplication_radix(result, X, Y, M, m, alpha, beta);

        // Print result
        gmp_printf("Result: %Zd\n", result);
        exit(0);
        // Cleanup
        mpz_clears(X, Y, M, result, NULL);



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
