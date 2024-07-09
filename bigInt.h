#include <stddef.h>
#include <stdint.h>

#ifndef CALCFIBONACCI_BIGINT_H
#define CALCFIBONACCI_BIGINT_H

struct bigInt {
    size_t size_32;
    uint32_t *ptr_32;
};
#endif //CALCFIBONACCI_BIGINT_H
