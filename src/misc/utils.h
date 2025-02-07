#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

uint64_t avx_sum(uint64_t *arr, const size_t len, int *carry_count);


#endif //UTILS_H
