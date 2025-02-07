#ifndef MOD_H
#define MOD_H

#include <gmp.h>
#include <stdint.h>

/// @brief Computing Remainder with respect to (2^64-1) and (2^64-2)
/// @param Number Input Number
/// @param sum64min1 Remainder w.r.t (2^64-1)
/// @param sum64min2 Remainder w.r.t (2^64-2)
void limb_sum_mod64m1_64m2(const mpz_t Number, uint64_t *sum64min1, uint64_t *sum64min2);

const static __uint128_t mod_64minus1 = ((__uint128_t) 1 << 64) - 1;
const static __uint128_t mod_64minus2 = ((__uint128_t) 1 << 64) - 2;


#endif
