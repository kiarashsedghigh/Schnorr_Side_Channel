#include "utils/mod.h"
#include <stdlib.h>
#include <stdio.h>

/* Global variables */
mpz_t r, mu, temp, mask, shift_amount;
mpz_t *Y_words;
uint64_t r_sum64min1, r_sum64min2;
uint64_t m_sum64min1, m_sum64min2;
mpz_t Z_curr, q, Z;
long d, N;
int Y_last_limb_size;

#define Y_WORDS_COUNT 256
#define BM_DEFEND_COMPUTE_SUM


void BM_encode_inputs(const mpz_t X, const mpz_t Y, mpz_t X_encoded, mpz_t Y_encoded, mpz_t M, const int shift) {
    /*
     * Encoding X =>    X = X + M
     * Encoding Y =>    Y = Y << shift
     */
    mpz_set(X_encoded, X);
    mpz_set(Y_encoded, Y);
    mpz_add(X_encoded, X_encoded, M);
    mpz_mul_2exp(Y_encoded, Y_encoded, shift);
}

void setup_BM_radix(const mpz_t M, const unsigned long m, const unsigned long alpha) {
    mpz_inits(r, mu, temp, mask, shift_amount, NULL);

    /* Precompute r = 2^m and mu = floor(2^(N+alpha) / M) */
    mpz_ui_pow_ui(r, 2, m); // r = 2^m
    N = mpz_sizeinbase(M, 2); // N = bit length of M
    d = (N + m - 1) / m; // ceiling(N / m)
    mpz_ui_pow_ui(temp, 2, N + alpha); // temp = 2^(N+alpha)
    mpz_fdiv_q(mu, temp, M); // mu = floor(temp / M)

    /* Initialize Y_words */
    Y_words = malloc(Y_WORDS_COUNT * sizeof(mpz_t));
    for (size_t i = 0; i < Y_WORDS_COUNT; i++)
        mpz_init(Y_words[i]);

    /* For encoding purposes */
    mpz_ui_pow_ui(mask, 2, m); // mask = 2^m
    mpz_sub_ui(mask, mask, 1); // mask = 2^m - 1

    /* Computing remainder of M w.r.t 2^64-1 and 2^64-2 */
    limb_sum_mod64m1_64m2(M, &m_sum64min1, &m_sum64min2);

    mpz_inits(Z_curr, q, Z, NULL);
}

void setup_clear_BM_radix() {
    free(Y_words);
    mpz_clears(shift_amount, r, mu, temp, mask, NULL);
}

void BM_radix(mpz_t result, mpz_t partial_computation_result, const mpz_t X, const mpz_t Y, const mpz_t M,
              const unsigned long m, const unsigned long alpha, const long beta, const int recompute_count) {
    /* Break Y into d Y_words */
    for (size_t i = 0; i < d; i++) {
        mpz_fdiv_r_2exp(Y_words[i], Y, m * (i + 1));
        mpz_fdiv_q_2exp(Y_words[i], Y_words[i], m * i);
    }

    uint64_t Z_old_sum64min1, Z_old_sum64min2;
    uint64_t Z_curr_sum64min1, Z_curr_sum64min2;
    uint64_t X_sum64min1, X_sum64min2;
    uint64_t Y_sum64min1, Y_sum64min2;
    uint64_t q_sum64min1, q_sum64min2;

#ifdef BM_DEFEND_COMPUTE_SUM
    __uint128_t fault_flag1 = 0, fault_flag2 = 0;

    Y_last_limb_size = mpz_sizeinbase(Y_words[d - 1], 2);

    /* Compute here for efficiency */
    limb_sum_mod64m1_64m2(X, &X_sum64min1, &X_sum64min2);
#endif

    for (long i = d - 1; i >= 0; i--) {

        mpz_mul(Z_curr, Z, r); // Z(i) = Z(i+1) * r
        mpz_addmul(Z_curr, X, Y_words[i]); // Z(i) += X * Yi

#ifdef BM_DEFEND_COMPUTE_SUM
        if (i < d - recompute_count) {
            /* BM Countermeasure */
            limb_sum_mod64m1_64m2(Z_curr, &Z_curr_sum64min1, &Z_curr_sum64min2);

            Y_sum64min1 = Y_words[i]->_mp_d[0] % mod_64minus1;
            Y_sum64min2 = Y_words[i]->_mp_d[0] % mod_64minus2;

            const __uint128_t mult_zr1 = Z_old_sum64min1;
            const __uint128_t mult_xy1 = (__uint128_t) X_sum64min1 * (__uint128_t) Y_sum64min1;
            if ((__uint128_t) Z_curr_sum64min1 >= ((__uint128_t) mult_zr1 + mult_xy1))
                fault_flag1 = ((__uint128_t) Z_curr_sum64min1 - ((__uint128_t) mult_zr1 + mult_xy1)) % mod_64minus1;
            else {
                fault_flag1 = ((mult_zr1 + mult_xy1) % mod_64minus1 - (__uint128_t) Z_curr_sum64min1) % mod_64minus1;
                if (fault_flag1 == 1)
                    fault_flag1 = 0;
            }

            const __uint128_t mult_zr2 = ((__uint128_t) Z_old_sum64min2) * 2;
            const __uint128_t mult_xy2 = ((__uint128_t) X_sum64min2 * (__uint128_t) Y_sum64min2);
            if ((__uint128_t) Z_curr_sum64min2 >= (mult_zr2 + mult_xy2))
                fault_flag2 = ((__uint128_t) Z_curr_sum64min2 - (mult_zr2 + mult_xy2) % mod_64minus2) % mod_64minus2;
            else
                fault_flag2 = ((mult_zr2 + mult_xy2) % mod_64minus2 - (__uint128_t) Z_curr_sum64min2) % mod_64minus2;
            if (fault_flag2 != 0)
                fault_flag2 = 0;
        }
#endif

        /* q(i) = floor((Z(i) / 2^(N+β)) * (μ / 2^(α-β))) */
        mpz_fdiv_q_2exp(temp, Z_curr, N + beta); // temp = floor(Z(i) / 2^(N+β))
        mpz_mul(temp, temp, mu); // temp *= μ
        mpz_fdiv_q_2exp(q, temp, alpha - beta); // q = temp / 2^(α-β)
        mpz_submul(Z_curr, q, M);
        mpz_set(Z, Z_curr);

        if (i < d - recompute_count) {
            limb_sum_mod64m1_64m2(q, &q_sum64min1, &q_sum64min2);
            Z_old_sum64min1 = ((__uint128_t) Z_curr_sum64min1 - (__uint128_t) q_sum64min1 * (__uint128_t) m_sum64min1) %
                              mod_64minus1;
            Z_old_sum64min2 = ((__uint128_t) Z_curr_sum64min2 - ((__uint128_t) q_sum64min2 * (__uint128_t) m_sum64min2))
                              % mod_64minus2;
        }
        /* Compute mod values for the next steps that we need a value fo Z_old_sum64min[1/2] */
        if (i == d - recompute_count)
            limb_sum_mod64m1_64m2(Z, &Z_old_sum64min1, &Z_old_sum64min2);

#ifdef BM_DEFEND_COMPUTE_SUM
        if (d - i == recompute_count)
            mpz_set(partial_computation_result, Z);

        if (fault_flag1 != 0 || fault_flag2 != 0) {
            printf("Fault Found \n");
            exit(0);
        }
#endif
    }
    if (mpz_cmp(Z, M) >= 0)
        mpz_sub(Z, Z, M);

    mpz_set(result, Z);
}


void BM_radix_recomputation(mpz_t result_partial, const mpz_t X, const mpz_t Y, const mpz_t M,
                            const unsigned long m, const unsigned long alpha, const long beta,
                            const int recompute_count, const int Y_bit_length) {
    /* Calculate Y_words[l-1] = Y >> (k - msb) */
    mpz_div_2exp(Y_words[recompute_count - 1], Y, Y_bit_length - Y_last_limb_size);

    unsigned int j = 1;
    for (int i = recompute_count - 2; i >= 0; i--) {
        /* Calculate shift_amount = k - msb - (m * j) */
        mpz_set_ui(shift_amount, Y_bit_length - Y_last_limb_size - m * j);

        /* Right shift Y by shift_amount */
        mpz_div_2exp(temp, Y, mpz_get_ui(shift_amount));

        /* Mask the result: temp = temp & mask */
        mpz_and(Y_words[i], temp, mask);
        j++;
    }

    mpz_init(Z);
    for (ssize_t i = recompute_count - 1; i >= 0; i--) {
        mpz_mul(Z_curr, Z, r); // Z(i) = Z(i+1) * r
        mpz_addmul(Z_curr, X, Y_words[i]); // Z(i) += X * Yi
        mpz_fdiv_q_2exp(temp, Z_curr, N + beta); // temp = floor(Z(i) / 2^(N+β))
        mpz_mul(temp, temp, mu); // temp *= μ
        mpz_fdiv_q_2exp(q, temp, alpha - beta); // q = temp / 2^(α-β)
        mpz_submul(Z_curr, q, M);
        mpz_set(Z, Z_curr);
    }
    mpz_set(result_partial, Z);
}
