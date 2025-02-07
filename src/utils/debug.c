#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>


void print_128(const __uint128_t num) {
    const uint64_t high = num >> 64;
    const uint64_t low = num;
    printf("Decimal: %" PRIu64 "%019" PRIu64 "\n", high, low);
}
