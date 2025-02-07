#include <immintrin.h>
#include <stdint.h>
#include <stdio.h>



void print_m256(const __m256i m256) {
    uint64_t result[4];
    _mm256_storeu_si256((__m256i*)result, m256);

    // Print the values
    printf("Result:\n");
    for (int i = 0; i < 4; i++) {
        printf("result[%d] = %llu\n", i, result[i]);
    }
}

static inline __m256i mm256_cmpge_epu64(const __m256i a, const __m256i b) {
    const __m256i sign_flip = _mm256_set1_epi64x(0x8000000000000000);
    __m256i a_signed = _mm256_xor_si256(a, sign_flip);
    __m256i b_signed = _mm256_xor_si256(b, sign_flip);
    return _mm256_or_si256(_mm256_cmpgt_epi64(a_signed, b_signed), _mm256_cmpeq_epi64(a_signed, b_signed));
}


// static inline __m256i mm256_cmpge_epu64(__m256i a, __m256i b) {
//     const __m256i sign_flip = _mm256_set1_epi64x(0x8000000000000000);
//     a = _mm256_xor_si256(a, sign_flip);
//     b = _mm256_xor_si256(b, sign_flip);
//     return _mm256_cmpgt_epi64(a, b) | _mm256_cmpeq_epi64(a, b);
// }


static inline uint64_t count_carry_horizontal_sum(__m256i a) {
    int carry = 0;
    uint64_t numbers[4];
    for (int i = 0; i < 4; i++)
        numbers[i] = _mm256_extract_epi64(a, i);  // Extract the i-th 64-bit value

    // Left
    uint64_t sum_left = numbers[0] + numbers[1];

    // Right
    uint64_t sum_right = numbers[2] + numbers[3];

    if (sum_left < numbers[0])
        carry++;

    if (sum_right < numbers[2])
        carry++;

    uint64_t sum_total = sum_left + sum_right;
    if (sum_total < sum_left)
        carry++;
    return carry;
}

uint64_t avx_sum(uint64_t *arr, const size_t len, int *carry_count) {
    __m256i accumulator = _mm256_setzero_si256();
    *carry_count = 0;

    for (size_t i = 0; i + 4 <= len; i += 4) {
        // Load 4 unsigned 64-bit integers into an AVX2 register
        const __m256i data_batch = _mm256_loadu_si256((__m256i *)&arr[i]);

        // Add numbers with accumulator
        accumulator = _mm256_add_epi64(accumulator, data_batch);

        // Check for carry
        const __m256i difference_sums = mm256_cmpge_epu64(accumulator, data_batch);

        uint64_t overflow_result[4];
        _mm256_storeu_si256(overflow_result, difference_sums);
        for (int j = 0; j < 4; j++)
            *carry_count += (overflow_result[j] == 0x0);
    }
    *carry_count += count_carry_horizontal_sum(accumulator);

    uint64_t accumulator_arr[4];
    _mm256_storeu_si256((__m256i*)accumulator_arr, accumulator);
    uint64_t sum = 0;
    for (int j = 0; j < 4; j++)
        sum += accumulator_arr[j];
    return sum;
}


// int main() {
//
//
//
//     uint64_t arr[] = {5781633234971700525, 5886778865283468297,    18007604941570302945,   6635100435205195395,
//                 0,                         0,                0,     0};
//     size_t len = sizeof(arr) / sizeof(arr[0]);
//     int overflow_count=0;
//
//     uint64_t sum = avx_sum(arr, len, &overflow_count);
//     printf("Sum: %llu\n", sum);
//     printf("OF: %llu\n", overflow_count);
//     return 0;
// }
