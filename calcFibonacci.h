#ifndef CALCFIBONACCI_CALCFIBONACCI_H
#define CALCFIBONACCI_CALCFIBONACCI_H

#include "bigInt.h"

enum mulMode {
    Square,
    Normal
};

struct matrix {
    struct bigInt f_n_nxt;
    struct bigInt f_n_curr;
    struct bigInt f_n_prv;
};

struct bigInt addNumber(const struct bigInt add1, const struct bigInt add2, size_t shift);

struct bigInt subNumber(const struct bigInt minuend, const struct bigInt subtrahend);

struct bigInt multiplyNumber(const struct bigInt a, const struct bigInt b);

struct bigInt karatsubaMulNumber(struct bigInt mul1, struct bigInt mul2);

struct bigInt karatsubaSquareNumber(struct bigInt mul);

struct matrix multiplyMatrix(struct matrix a, struct matrix b, enum mulMode mode);

void clearMatrix(struct matrix m);

void fib(uint64_t n, size_t len, uint8_t buf[len]);

#endif //CALCFIBONACCI_CALCFIBONACCI_H