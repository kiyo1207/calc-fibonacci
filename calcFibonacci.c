#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <printf.h>
#include "calcFibonacci.h"


/**
 * Addition of add1 and left shifted add2(add1 + (add2 << shift)).
 * @param add1
 * @param add2
 * @return Result of addition
 */
struct bigInt addNumber(const struct bigInt add1, const struct bigInt add2, size_t shift) {

    size_t size = add1.size_32 >= add2.size_32 + shift ? add1.size_32 : add2.size_32 + shift;
    uint32_t *array;
    array = malloc(size * sizeof(uint32_t));
    if (array == NULL) {
        printf("Couldn't allocate memory : addNumber().\n");
        abort();
    }
    size_t index = 0;
    uint8_t carry = 0;

    //sum of the first shift-bits are fixed.
    while (index < shift) {
        array[index] = (index < add1.size_32) ? add1.ptr_32[index] : 0;
        index++;
    }

    while (index < add1.size_32 && index < add2.size_32 + shift) {
        uint64_t sum = (uint64_t) add1.ptr_32[index] + (uint64_t) add2.ptr_32[index - shift] + carry;
        array[index] = sum & 0xffffffffffffffff;
        carry = (sum >> 32);
        index++;
    }

    while (index < add1.size_32) {
        uint64_t sum = (uint64_t) add1.ptr_32[index] + carry;
        array[index] = sum & 0xffffffffffffffff;
        carry = (sum >> 32);
        index++;
    }

    while (index < add2.size_32 + shift) {
        uint64_t sum = (uint64_t) add2.ptr_32[index - shift] + carry;
        array[index] = sum & 0xffffffffffffffff;
        carry = (sum >> 32);
        index++;
    }

    if (carry) {
        array = realloc(array, (size + 1) * sizeof(uint32_t));
        array[size] = 1;
        return ((struct bigInt) {.size_32 = size + 1, .ptr_32 =  array});
    }
    return ((struct bigInt) {.size_32 = size, .ptr_32 =  array});
}


/**
 * Subtraction of two bigInts.
 * @param minuend
 * @param subtrahend
 * @return Result of subtraction
 */
struct bigInt subNumber(const struct bigInt minuend, const struct bigInt subtrahend) {

    size_t size = minuend.size_32;
    uint32_t *array = malloc(size * sizeof(uint32_t));
    if (array == NULL) {
        printf("Couldn't allocate memory : subNumber().\n");
        abort();
    }

    size_t index = 0;
    uint8_t carry = 1;

    // subtract by using 2's complement.
    while (index < minuend.size_32 && index < subtrahend.size_32) {
        uint64_t sum = (uint64_t) *(minuend.ptr_32 + index) + (uint64_t) ~*(subtrahend.ptr_32 + index) + carry;
        *(array + index) = sum;
        carry = (sum >> 32);
        index++;
    }

    while (index < minuend.size_32) {
        uint64_t sum = (uint64_t) *(minuend.ptr_32 + index) + carry + 0xffffffff;
        *(array + index) = sum;
        carry = (sum >> 32);
        index++;
    }

    size_t modified_size = size;
    while (modified_size > 1
           && array[modified_size - 1] == 0)
        modified_size--;

    if (modified_size < size) {
        array = realloc(array, modified_size * sizeof(uint32_t));
        return ((struct bigInt) {.size_32 = modified_size, .ptr_32 =  array});
    }
    return ((struct bigInt) {.size_32 = size, .ptr_32 =  array});

}


/**
 * Multiply two bigInt.
 * @param a
 * @param b
 * @return Result of subtraction
 */
struct bigInt multiplyNumber(struct bigInt a, struct bigInt b) {

    // Calculate maximum length
    size_t maxLength = a.size_32 + b.size_32;

    uint64_t tmpProduct[maxLength];
    for (size_t i = 0; i < maxLength; i++) {
        tmpProduct[i] = 0;
    }

    // For memory access efficiency
    if (a.size_32 > b.size_32) {
        struct bigInt tmp = a;
        a = b;
        b = tmp;
    }

    for (size_t i = 0; i < a.size_32; i++) {
        for (size_t j = 0; j < b.size_32; j++) {
            tmpProduct[i + j] += (uint64_t) a.ptr_32[i] * (uint64_t) b.ptr_32[j];
            size_t k = i + j;
            while (tmpProduct[k] > UINT32_MAX) {
                tmpProduct[k + 1] += tmpProduct[k] >> 32;
                tmpProduct[k] &= 0xFFFFFFFF;
                k++;
            }
        }
    }

    // Shorten leading zeros
    size_t shortenedLength = maxLength;
    while (shortenedLength > 1 && tmpProduct[shortenedLength - 1] == 0) {
        shortenedLength--;
    }

    struct bigInt ret;
    ret.size_32 = shortenedLength;
    ret.ptr_32 = malloc(shortenedLength * sizeof(uint32_t));
    if (ret.ptr_32 == NULL) {
        printf("Couldn't allocate memory : multiplyNumber().\n");
        abort();
    }

    for (size_t i = 0; i < ret.size_32; i++) {
        ret.ptr_32[i] = (uint32_t) tmpProduct[i];
    }

    return ret;
}


/**
 * Square bigInt.
 * @param a
 * @return Result of subtraction
 */
struct bigInt squareNumber(struct bigInt a) {

    // Calculate maximum length
    size_t maxLength = a.size_32 * 2;

    uint64_t tmpProduct[maxLength];
    for (size_t i = 0; i < maxLength; i++) {
        tmpProduct[i] = 0;
    }

    for (int i = 0; i < a.size_32; i++) {
        for (int j = i; j < a.size_32; j++) {
            uint64_t tmp = (uint64_t) a.ptr_32[i] * (uint64_t) a.ptr_32[j];
            for (int k = 0; k < 2 - (i == j); k++) {
                tmpProduct[i + j] += tmp;
                size_t l = i + j;
                while (tmpProduct[l] > UINT32_MAX) {
                    tmpProduct[l + 1] += tmpProduct[l] >> 32;
                    tmpProduct[l] &= 0xFFFFFFFF;
                    l++;
                }
            }
        }
    }

    // Shorten leading zeros
    size_t shortenedLength = maxLength;
    while (shortenedLength > 1 && tmpProduct[shortenedLength - 1] == 0) {
        shortenedLength--;
    }

    struct bigInt ret;
    ret.size_32 = shortenedLength;
    ret.ptr_32 = malloc(shortenedLength * sizeof(uint32_t));
    if (ret.ptr_32 == NULL) {
        printf("Couldn't allocate memory : squareNumber().\n");
        abort();
    }

    for (size_t i = 0; i < ret.size_32; i++) {
        ret.ptr_32[i] = (uint32_t) tmpProduct[i];
    }

    return ret;
}

/**
 * Multiply two bigInts with Karatsuba algorithm.
 * Apply the Karatsuba algorithm recursively until the sum of the length of mul1 and mul2 is less than 50.
 * @param mul1
 * @param mul2
 * @return Result of multiplication
 */
struct bigInt karatsubaMulNumber(struct bigInt mul1, struct bigInt mul2) {


    if (mul1.size_32 + mul2.size_32 <= 50) {
        return multiplyNumber(mul1, mul2);
    }

    if (mul1.size_32 < mul2.size_32) {
        struct bigInt tmp = mul1;
        mul1 = mul2;
        mul2 = tmp;
    }

    size_t shift = mul1.size_32 / 2;

    struct bigInt mul1_lower;
    mul1_lower.size_32 = shift;
    mul1_lower.ptr_32 = mul1.ptr_32;


    struct bigInt mul1_upper;
    mul1_upper.size_32 = mul1.size_32 - shift;
    mul1_upper.ptr_32 = mul1.ptr_32 + shift;

    if (mul2.size_32 <= mul1_lower.size_32) {
        struct bigInt u1_mul_mul2 = karatsubaMulNumber(mul1_upper, mul2);
        struct bigInt l1_mul_mul2 = karatsubaMulNumber(mul1_lower, mul2);
        struct bigInt u1_mul_mul2_add_l1_mul_mul2 = addNumber(l1_mul_mul2, u1_mul_mul2, shift);
        free(u1_mul_mul2.ptr_32);
        free(l1_mul_mul2.ptr_32);
        return u1_mul_mul2_add_l1_mul_mul2;
    }

    struct bigInt mul2_lower;
    mul2_lower.size_32 = shift;
    mul2_lower.ptr_32 = mul2.ptr_32;
    struct bigInt mul2_upper;
    mul2_upper.size_32 = mul2.size_32 - shift;
    mul2_upper.ptr_32 = mul2.ptr_32 + shift;

    struct bigInt l1_mul_l2 = karatsubaMulNumber(mul1_lower, mul2_lower);
    struct bigInt u1_mul_u2 = karatsubaMulNumber(mul1_upper, mul2_upper);

    struct bigInt u1_add_l1 = addNumber(mul1_upper, mul1_lower, 0);
    struct bigInt u2_add_l2 = addNumber(mul2_upper, mul2_lower, 0);

    struct bigInt u1_add_l1_mul_u2_add_l2 = karatsubaMulNumber(u1_add_l1, u2_add_l2);

    struct bigInt subtrahend = addNumber(u1_mul_u2, l1_mul_l2, 0);
    struct bigInt u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2 = subNumber(u1_add_l1_mul_u2_add_l2, subtrahend);

    struct bigInt tmp = addNumber(l1_mul_l2, u1_mul_u2, shift * 2);
    struct bigInt sum = addNumber(tmp, u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2,
                                  shift);

    free(tmp.ptr_32);
    free(u1_add_l1_mul_u2_add_l2.ptr_32);
    free(u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2.ptr_32);
    free(subtrahend.ptr_32);
    free(u1_add_l1.ptr_32);
    free(u2_add_l2.ptr_32);
    free(l1_mul_l2.ptr_32);
    free(u1_mul_u2.ptr_32);

    return sum;

}

/**
 * square bigInt with Karatsuba algorithm.
 * Apply the Karatsuba algorithm recursively until the length of mul is less than 50.
 * @param mul1
 * @param mul2
 * @return Result of square
 */
struct bigInt karatsubaSquareNumber(struct bigInt mul) {

    if (mul.size_32 <= 25) {
        return squareNumber(mul);
    }

    size_t shift = mul.size_32 / 2;

    struct bigInt mul_lower;
    mul_lower.size_32 = shift;
    mul_lower.ptr_32 = mul.ptr_32;

    struct bigInt mul_upper;
    mul_upper.size_32 = mul.size_32 - shift;
    mul_upper.ptr_32 = mul.ptr_32 + shift;

    struct bigInt squareL = karatsubaSquareNumber(mul_lower);
    struct bigInt squareU = karatsubaSquareNumber(mul_upper);
    struct bigInt addUL = addNumber(mul_upper, mul_lower, 0);
    struct bigInt square_addUL = karatsubaSquareNumber(addUL);
    struct bigInt add_squareL_squareU = addNumber(squareL, squareU, 0);
    struct bigInt u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2 = subNumber(square_addUL, add_squareL_squareU);

    struct bigInt tmp = addNumber(squareL, squareU, shift * 2);
    struct bigInt sum = addNumber(tmp, u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2,
                                  shift);

    free(tmp.ptr_32);
    free(square_addUL.ptr_32);
    free(u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2.ptr_32);
    free(add_squareL_squareU.ptr_32);
    free(addUL.ptr_32);
    free(squareU.ptr_32);
    free(squareL.ptr_32);

    return sum;
}


/**
 * Multiply two 2*2 matrices. Each matrix has (n+1)th, nth,(n-1)th elements.
 * Multiply　mode can be chosen, the same matrix multiplication (square) or normal multiplication.
 *
 * @param a matrix(f(n+1), f(n),f(n-1))
 * @param b matrix(f(n+1), f(n),f(n-1))
 * @param mode when Normal is set then normal matrix multiplication, when Square then square.
 * @return product of 2 Matrices.
 */
struct matrix multiplyMatrix(struct matrix a, struct matrix b, enum mulMode mode) {

    struct bigInt nxt_mul_nxt;
    struct bigInt n_mul_n;
    struct bigInt pre_mul_pre;
    struct matrix m;

    if (mode == Normal) {
        // normal matrix multiplication
        nxt_mul_nxt = karatsubaMulNumber(a.f_n_nxt, b.f_n_nxt);
        n_mul_n = karatsubaMulNumber(a.f_n_curr, b.f_n_curr);
        pre_mul_pre = karatsubaMulNumber(a.f_n_prv, b.f_n_prv);
        m.f_n_nxt = addNumber(nxt_mul_nxt, n_mul_n, 0);
        m.f_n_prv = addNumber(n_mul_n, pre_mul_pre, 0);
        m.f_n_curr = subNumber(m.f_n_nxt, m.f_n_prv);
    } else {
        // square matrix multiplication
        nxt_mul_nxt = karatsubaSquareNumber(a.f_n_nxt);
        n_mul_n = karatsubaSquareNumber(a.f_n_curr);
        pre_mul_pre = karatsubaSquareNumber(a.f_n_prv);
        m.f_n_nxt = addNumber(nxt_mul_nxt, n_mul_n, 0);
        m.f_n_prv = addNumber(n_mul_n, pre_mul_pre, 0);
        m.f_n_curr = subNumber(nxt_mul_nxt, pre_mul_pre);
    }

    free(n_mul_n.ptr_32);
    free(pre_mul_pre.ptr_32);
    free(nxt_mul_nxt.ptr_32);

    return m;
}


void clearMatrix(struct matrix m) {
    free(m.f_n_prv.ptr_32);
    free(m.f_n_curr.ptr_32);
    free(m.f_n_nxt.ptr_32);
}


void fib(uint64_t n, size_t len, uint8_t buf[len]) {

    uint32_t *squareNxt = malloc(sizeof(uint32_t));
    uint32_t *squareCurr = malloc(sizeof(uint32_t));
    uint32_t *squarePrv = calloc(1, sizeof(uint32_t));
    uint32_t *resultNxt = malloc(sizeof(uint32_t));
    uint32_t *resultCurr = calloc(1, sizeof(uint32_t));
    uint32_t *resultPrv = malloc(sizeof(uint32_t));

    if (squareNxt == NULL || squareCurr == NULL || squarePrv == NULL || resultNxt == NULL || resultCurr == NULL ||
        resultPrv == NULL) {
        printf("Couldn't allocate memory : fib().\n");
        abort();
    }
    *(squareNxt) = 1;
    *(squareCurr) = 1;
    *(resultNxt) = 1;
    *(resultPrv) = 1;


    struct matrix squareMatrix;
    squareMatrix.f_n_prv.ptr_32 = squarePrv;
    squareMatrix.f_n_prv.size_32 = 1;
    squareMatrix.f_n_curr.ptr_32 = squareCurr;
    squareMatrix.f_n_curr.size_32 = 1;
    squareMatrix.f_n_nxt.ptr_32 = squareNxt;
    squareMatrix.f_n_nxt.size_32 = 1;

    // identity matrix is the basic.
    struct matrix resultMatrix;
    resultMatrix.f_n_prv.ptr_32 = resultPrv;
    resultMatrix.f_n_prv.size_32 = 1;
    resultMatrix.f_n_curr.ptr_32 = resultCurr;
    resultMatrix.f_n_curr.size_32 = 1;
    resultMatrix.f_n_nxt.ptr_32 = resultNxt;
    resultMatrix.f_n_nxt.size_32 = 1;

    //calculate fib(size_32)
    while (n > 0) {
        //multiply res with square
        if (n & 0x1) {
            struct matrix newResultMatrix = multiplyMatrix(resultMatrix, squareMatrix, Normal);
            clearMatrix(resultMatrix);
            resultMatrix = newResultMatrix;
        }

        //Termination condition
        n = n >> 1;
        if (!n) break;
        struct matrix newSquareMatrix = multiplyMatrix(squareMatrix, squareMatrix, Square);
        clearMatrix(squareMatrix);
        squareMatrix = newSquareMatrix;
    }

    //Write result to the buffer
    size_t i = 0;
    for (; i < resultMatrix.f_n_curr.size_32 * 4; i += 4) {
        uint32_t val = resultMatrix.f_n_curr.ptr_32[i / 4];
        buf[i + 0] = (val >> 0) & 0xff;
        buf[i + 1] = (val >> 8) & 0xff;
        buf[i + 2] = (val >> 16) & 0xff;
        buf[i + 3] = (val >> 24) & 0xff;
    }

    //padding zero
    while (i < len) {
        buf[i] = 0;
        i++;
    }
    clearMatrix(resultMatrix);
    clearMatrix(squareMatrix);
}
