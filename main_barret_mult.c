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

        mpz_t X, Y, M, result;
        mpz_inits(X, Y, M, result, NULL);

        // Initialize values for X, Y, M
        mpz_set_str(X, "20214201316202783299242554811138106343922050269439575514811953896596770920561145830633211474257353477419096402680177431711860578647926876212828203432964325144562947296861436637069154268326440130616295489416570794279766291397798123520937149287459130764295397249712099718426268653072941374565881877599757060485180302111477255959027025317470927602078073284308908529611367116151135474314589582021784465957342854085667419164935035077096607128269228716598056975187794345076923710655402498715478026344484845498671904503044434165072041795123288486395840986646767256168515806172446481341106247703034104109062349182746850571900", 10);
        mpz_set_str(Y, "16832320556202393978997812802647765108318837276732427572851460854482540460463531843906520611803832294408098040718124280854840056228220581493908818324241124188155145958367050420768910508436935872893861725146421501055408547373761685242498971431237054627158088945749870468726167188286543345911924152914141659317586936725441711836455784582934591821100342178381596214498268839737163666111726461227024540370351399002476556840267673870968380701370011161193466493138715720715562410586014527627171277676285268375632564131544090851431264976552774929982791110994661428824452409378975591449782014864471863746609062880608030030940", 10);
        mpz_set_str(M, "86832320556202393978997812802647765108318837276732427572851460854482540460463531843906520611803832294408098040718124280854840056228220581493908818324241124188155145958367050420768910508436935872893861725146421501055408547373761685242498971431237054627158088945749870468726167188286543345911924152914141659317586936725441711836455784582934591821100342178381596214498268839737163666111726461227024540370351399002476556840267673870968380701370011161193466493138715720715562410586014527627171277676285268375632564131544090851431264976552774929982791110994661428824452409378975591449782014864471863746609062880608030030940", 10);



        // Perform Barrett Modular Multiplication
        unsigned long m = 32; // word size
        long alpha = m+3, beta = -2; // parameters
        barret_multiplication_radix(result, X, Y, M, m, alpha, beta);

        // Print result
        gmp_printf("Result_our: %Zd\n", result);

        mpz_mul(X,X,Y);
        mpz_mod(X,X,M);


        gmp_printf("Result GMP: %Zd\n", X);

        mpz_sub(X,X,result);
        gmp_printf("Result EQ: %Zd\n", X);



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
