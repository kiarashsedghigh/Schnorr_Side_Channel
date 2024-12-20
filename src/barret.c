#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "size.h"
#include <stdint.h>

//
//
// __uint128_t mod_64plus1 = ((__uint128_t)1 << 64) + 1; // Use 128-bit integer
__uint128_t mod_64minus1 = ((__uint128_t)1 << 64) - 1; // Use 128-bit integer
__uint128_t mod_64minus2 = ((__uint128_t)1 << 64) - 2; // Use 128-bit integer
//
// __uint128_t two_256 = ((__uint128_t)1 << 256);
// __uint128_t mu = ((__uint128_t)1 << 256) / mod_64plus1;
//
// void print_uint128(__uint128_t value) {
//     uint64_t high = (uint64_t)(value >> 64); // Extract the high 64 bits
//     uint64_t low = (uint64_t)value;         // Extract the low 64 bits
//
//     if (high > 0) {
//         // If the high part is non-zero, print it with leading zeros for the low part
//         printf("High64:%llu      Low64:%llu\n", high, low);
//     } else {
//         // If the high part is zero, print only the low part
//         printf("High64: 0        Low64: %llu\n", low);
//     }
// }
//
// uint64_t sum_eff_mod64min1(mpz_t N) {
//     const mp_size_t limb_cnt = mpz_size(N);  // Number of limbs in the mpz_t number
//     const mp_limb_t *limbs = mpz_limbs_read(N);  // Get the limbs of the number
//
//     __uint128_t temp_sum = 0;
//
//     for (int i = 0; i < limb_cnt; i++) {
//         temp_sum += limbs[i];
//         if (temp_sum >= mod_64minus1) {
//             temp_sum -= mod_64minus1;
//         }
//     }
//     return temp_sum;
// }
//
//
//
// __uint128_t N64MIN2;// = ((__uint128_t)1 << 64) - 2;
// __uint128_t N64MIN1;// = ((__uint128_t)1 << 64) - 2;
//
// __uint128_t B64;// = ((__uint128_t)1 << 64);
// __uint128_t MOD65MIN1;// = ((__uint128_t)1 << 65) - 1;
// __uint128_t POW65;// = ((__uint128_t)1 << 65);
// __uint128_t R;
// __uint128_t R2;
//
//
// void init_r() {
//     N64MIN2 = ((__uint128_t)1 << 64) - 2;
//     B64 = ((__uint128_t)1 << 64);
//     MOD65MIN1 = ((__uint128_t)1 << 65) - 1;
//     POW65 = ((__uint128_t)1 << 65);
//     R = ((__uint128_t)1 << 64) | 2;
//
//
//     N64MIN1 = ((__uint128_t)1 << 64) - 1;
//     R2 = ((__uint128_t)1 << 64) | 1;
//
// }
//
// __uint128_t barret_64min2_efficient_version(__uint128_t u) {
//     __uint128_t q_hat = (((u >> 63) * R) >> 65);
//     __uint128_t r1 = u & MOD65MIN1;
//     __uint128_t r2 = (q_hat * N64MIN2) & MOD65MIN1;
//
//     __uint128_t r;
//
//     if (r1 < r2) {
//         r = POW65 + r1 - r2;
//     }else {
//         r= r1 - r2;
//     }
//     while (r >= N64MIN2) {
//         r -= N64MIN2;
//     }
//     return r;
// }
//
//
// __uint128_t barret_64min1_efficient_version(__uint128_t u) {
//     printf("Value of u\n");
//     print_uint128(u);
//     __uint128_t q_hat = (((u >> 63) * R2) >> 65);
//
//     __uint128_t r1 = u & MOD65MIN1;
//     __uint128_t r2 = (q_hat * N64MIN1) & MOD65MIN1;
//
//
//     __uint128_t r;
//
//     if (r1 < r2) {
//         r = POW65 + r1 - r2;
//     }else {
//         r= r1 - r2;
//     }
//
//     while (r >= N64MIN1) {
//         r -= N64MIN1;
//     }
//
//     // Step 9: Return the result
//     return r;
// }
//
//
void sum_eff_mod64twomodes(mpz_t N, uint64_t* sum64min1, uint64_t* sum64min2) {
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
//
//
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

//
// __uint128_t barrett_reduction_128_64plus1(__uint128_t x) {
//     __uint128_t mu, q, r;
//
//     // Compute q = (x * mu) >> 128
//     q = (x * mu) >> 128;
//
//     // Compute r = x - q * N
//     r = x - (q * mod_64plus1);
//
//     // If r >= N, subtract N
//     if (r >= mod_64plus1) {
//         r -= mod_64plus1;
//     }
//     return r;
// }
//
// uint64_t m = (1ULL << 64) - 1;
//
// // Barrett reduction for modulus 2^64 - 1
// uint64_t barrett_reduction_128_64minus1(__uint128_t a) {
//     // m = 2^64 - 1, so we need to reduce a modulo (2^64 - 1)
//
//     // Compute q = floor(a / 2^64)
//     uint64_t q = (uint64_t)(a >> 64);
//
//     // Compute r = a - q * (2^64 - 1)
//     uint64_t r = (uint64_t)(a & m) - q;
//
//     // If r >= m, subtract m to get the final result
//     if (r >= m) {
//         r -= m;
//     }
//
//     return r;
// }
//
//
// // void barrett_reduction_64plus1(mpz_t N) {
//
//
//     // mpz_t u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2;
//     //
//     // // Initialize temporary variables
//     // mpz_inits(u_high, q_hat, r1, r2, temp, N_mult, b_n_plus_2, NULL);
//     //
//     // // Step 1: Compute q_hat = floor(floor(u / b^(n-1)) * R / b^(n+1))
//     // mpz_tdiv_q_2exp(u_high, u, (n - 1) * WORD_SIZE); // u_high = u >> (n-1)*word_size
//     // mpz_mul(q_hat, u_high, R);                       // q_hat = u_high * R
//     // mpz_tdiv_q_2exp(q_hat, q_hat, (n + 1) * WORD_SIZE); // q_hat = q_hat >> (n+1)*word_size
//     //
//     // // Step 2: Compute r = u mod b^(n+1) - (q_hat * N) mod b^(n+1)
//     // // Instead of using mpz_mod, use bitwise AND
//     // mpz_sub_ui(b_n_plus_2, b_n_plus_1, 1); // b_n_plus_1 - 1
//     //
//     // // r1 = u & (b_n_plus_1 - 1) (equivalent to u % b^(n+1))
//     // mpz_and(r1, u, b_n_plus_2);
//     //
//     // // r2 = (q_hat * N) & (b_n_plus_1 - 1) (equivalent to (q_hat * N) % b^(n+1))
//     // mpz_mul_ui(N_mult, q_hat, mod_64plus1);  // N_mult = q_hat * N
//     // mpz_and(r2, N_mult, b_n_plus_2);  // r2 = (q_hat * N) & (b_n_plus_1 - 1)
//     //
//     // // r = r1 - r2
//     // mpz_sub(r, r1, r2);
//     //
//     // // Step 3: If r < 0, add b^(n+1)
//     // if (mpz_sgn(r) < 0) {
//     //     mpz_add(r, r, b_n_plus_2);
//     // }
//     //
//     // // Step 4: While r >= N, subtract N
//     // while (mpz_cmp_ui(r, mod_64plus1) >= 0) {
//     //     mpz_sub_ui(r, r, mod_64plus1);
//     // }
//     //
//     // mpz_clears(u_high, q_hat, r1, r2, temp, N_mult, NULL);
// // }
//
