#include <gmp.h>
#include <stdint.h>
#include "mod.h"


void limb_sum_mod64m1_64m2(const mpz_t Number, uint64_t *sum64min1, uint64_t *sum64min2) {
    const mp_size_t limb_cnt = Number->_mp_size;
    const mp_limb_t *limbs = Number->_mp_d;
    __uint128_t temp_sum_64min1 = 0;
    __uint128_t temp_sum_64min2 = 0;

    for (int i = 0; i < limb_cnt; i++) {
        const __uint128_t limb = limbs[i];
        /* 2^64 - 1 */
        temp_sum_64min1 += limb;

        /* 2^64 - 2 */
        temp_sum_64min2 += limb << i;
    }
    // ((uint64_t *)&temp_sum_64min1)[1] + (temp_sum_64min1 != mod_64minus1) * (uint64_t)temp_sum_64min1;
    *sum64min1 = temp_sum_64min1 % mod_64minus1;
    *sum64min2 = temp_sum_64min2 % mod_64minus2;
}
