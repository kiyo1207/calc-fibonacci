//
// Created by Kyosuke Nakanishi on 08.07.24.
//

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <printf.h>
#include "calcFibonacci.h"




//////////
/**
 * addition of two bigInts.
 * @param add1
 * @param add2
 * @return Result of addition
 */
struct bigInt addNumber(const struct bigInt add1, const struct bigInt add2, size_t shift) {

    size_t size = add1.size_32 >= add2.size_32 + shift ? add1.size_32 : add2.size_32 + shift;
    uint32_t *array;
    array = malloc(size * sizeof(uint32_t));
    size_t index = 0;
    uint8_t carry = 0;

    //sum of the first shift-bits are fixed.
    while (index < shift) {
        array[index] = (index < add1.size_32) ? add1.ptr_32[index] : 0;
        index++;
    }

    while (index < add1.size_32 && index < add2.size_32 + shift) {
        uint64_t sum = ((uint64_t) add1.ptr_32[index]) + ((uint64_t) add2.ptr_32[index - shift]) + carry;
        array[index] = sum & 0xffffffffffffffff;
        carry = (sum >> 32);
        index++;
    }

    while (index < add1.size_32) {
        uint64_t sum = ((uint64_t) add1.ptr_32[index]) + carry;
        array[index] = sum & 0xffffffffffffffff;
        carry = (sum >> 32);
        index++;
    }

    while (index < add2.size_32 + shift) {
        uint64_t sum = ((uint64_t) add2.ptr_32[index - shift]) + carry;
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
 * subtraction of two bigInts.
 * @param minuend
 * @param subtrahend
 * @return Result of subtraction
 */
struct bigInt subNumber(const struct bigInt minuend, const struct bigInt subtrahend) {

    size_t size = minuend.size_32;
    uint32_t *array = malloc(size * sizeof(uint32_t));

    if (array == 0) {
        printf("Could not reserve memory for subtraction.\n");
        abort();
    }

    size_t index = 0;
    uint8_t carry = 1;

    // subtract by using 2's complement.
    while (index < minuend.size_32 && index < subtrahend.size_32) {
        uint64_t sum = ((uint64_t) *(minuend.ptr_32 + index)) + ((uint64_t) ~*(subtrahend.ptr_32 + index)) + carry;
        *(array + index) = sum;
        carry = (sum >> 32);
        index++;
    }

    while (index < minuend.size_32) {
        uint64_t sum = ((uint64_t) *(minuend.ptr_32 + index)) + carry + 0xffffffff;
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


struct bigInt multiplyNumber(const struct bigInt a, const struct bigInt b) {
    // Calculate maximum length
    size_t maxLength = a.size_32 + b.size_32;

    // Zero the intermediate result
    uint64_t tmpProduct[maxLength];
    for (size_t i = 0; i < maxLength; i++) {
        tmpProduct[i] = 0;
    }

    // Perform written multiplication
    for (size_t i = 0; i < a.size_32; i++) {
        for (size_t j = 0; j < b.size_32; j++) {
            tmpProduct[i + j] += (uint64_t) a.ptr_32[i] * (uint64_t) b.ptr_32[j];

            // Handle overflow
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

    // Copy result to return array
    struct bigInt ret;
    ret.size_32 = shortenedLength;
    ret.ptr_32 = malloc(shortenedLength * sizeof(uint32_t));
    if (ret.ptr_32 == NULL) {
        printf("Couldn't allocate memory in mul.\n");
        abort();
    }

    for (size_t i = 0; i < ret.size_32; i++) {
        ret.ptr_32[i] = (uint32_t) tmpProduct[i];
    }

    return ret;
}

struct bigInt karatubaMulNumber(struct bigInt mul1, struct bigInt mul2) {

    if (mul1.size_32 + mul2.size_32 <= 50) {
        return multiplyNumber(mul1, mul2);
    }

    long mul1_index_2nd = (mul1.size_32 / 2);
    long mul2_index_2nd = (mul2.size_32 / 2);

    if (mul1_index_2nd < mul2_index_2nd) {
        struct bigInt tmp = mul1;
        mul1 = mul2;
        mul2 = tmp;
        long tmp_index = mul1_index_2nd;
        mul1_index_2nd = mul2_index_2nd;
        mul2_index_2nd = tmp_index;
    }

    struct bigInt mul1_lower;
    mul1_lower.size_32 = mul1.size_32 / 2;
    mul1_lower.ptr_32 = mul1.ptr_32;

    struct bigInt mul1_upper;
    mul1_upper.size_32 = mul1.size_32 - mul1_index_2nd;
    mul1_upper.ptr_32 = &mul1.ptr_32[mul1_index_2nd];
    long mul1_upper_shift = mul1_index_2nd;

    struct bigInt mul2_lower;
    mul2_lower.size_32 = mul2.size_32 / 2;
    mul2_lower.ptr_32 = mul2.ptr_32;

    struct bigInt mul2_upper;
    mul2_upper.size_32 = mul2.size_32 - mul2_index_2nd;
    mul2_upper.ptr_32 = &mul2.ptr_32[mul2_index_2nd];
    long mul2_upper_shift = mul2_index_2nd;


    struct bigInt l1_mul_l2 = karatubaMulNumber(mul1_lower, mul2_lower);
    struct bigInt u1_mul_u2 = karatubaMulNumber(mul1_upper, mul2_upper);

    long shift1 = 0;
    long shift2 = 0;
    long shift3 = mul2_index_2nd;
    if (mul1_index_2nd != mul2_index_2nd) {
        shift1 = mul1_index_2nd - mul2_index_2nd;
        shift2 = mul1_index_2nd - mul2_index_2nd;
    }

    struct bigInt u1_add_l1 = addNumber(mul1_upper, mul1_lower, 0);
    struct bigInt u2_add_l2 = addNumber(mul2_upper, mul2_lower, shift1);
    struct bigInt u1_add_l1_mul_u2_add_l2 = karatubaMulNumber(u1_add_l1, u2_add_l2);
    struct bigInt subtrahend = addNumber(u1_mul_u2, l1_mul_l2, shift2);
    struct bigInt u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2 = subNumber(u1_add_l1_mul_u2_add_l2, subtrahend);

    struct bigInt sum1 = l1_mul_l2;
    struct bigInt sum1_sum2 = addNumber(sum1, u1_mul_u2, mul1_upper_shift + mul2_upper_shift);
    struct bigInt sum1_sum2_sum3 = addNumber(sum1_sum2, u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2, shift3);

    free(sum1_sum2.ptr_32);
    free(u1_add_l1_mul_u2_add_l2.ptr_32);
    free(u1_add_l1_mul_u2_add_l2_sub_l1_mul_l2_u1_mul_u2.ptr_32);
    free(subtrahend.ptr_32);
    free(u1_add_l1.ptr_32);
    free(u2_add_l2.ptr_32);
    free(u1_mul_u2.ptr_32);

    return sum1_sum2_sum3;
}

/**
 * multiply supermassive two same Number(bigInt) with Karatsuba algorithm.
 * continue to apply the Karatsuba algorithm recursively until the length of either mul1 or mul2 is 1.
 * @param mul1
 * @param mul2
 * @return answer of multiplication
 */
struct bigInt karatsubaMulSameNumber(struct bigInt mul) {

    if (mul.size_32 <= 25) {
        return multiplyNumber(mul, mul);
    }

    long mul_index = (mul.size_32 / 2);

    // lower parts of mul
    struct bigInt mul_lower;
    mul_lower.size_32 = mul.size_32 / 2;
    mul_lower.ptr_32 = mul.ptr_32;

    // upper parts of mul
    struct bigInt mul_upper;
    mul_upper.size_32 = mul.size_32 - mul_index;
    mul_upper.ptr_32 = &mul.ptr_32[mul_index];
    long mul_upper_shift = mul_index;

    struct bigInt l_mul_l = karatsubaMulSameNumber(mul_lower);
    struct bigInt u_mul_u = karatsubaMulSameNumber(mul_upper);

    long shift3 = mul_index;
    struct bigInt u_add_l = addNumber(mul_upper, mul_lower, 0);
    struct bigInt u_add_l_mul_u_add_l = karatsubaMulSameNumber(u_add_l);
    struct bigInt subtrahend = addNumber(u_mul_u, l_mul_l, 0);
    struct bigInt u_add_l_mul_u_add_l_sub_l_mul_l_u_mul_u = subNumber(u_add_l_mul_u_add_l, subtrahend);

    struct bigInt sum1 = l_mul_l;
//    struct bigInt sum2 = shift(u_mul_u, mul_upper_shift * 2);
//    struct bigInt sum3 = shift(u_add_l_mul_u_add_l_sub_l_mul_l_u_mul_u, shift3);
    struct bigInt sum1_sum2 = addNumber(sum1, u_mul_u, mul_upper_shift * 2);
    struct bigInt sum1_sum2_sum3 = addNumber(sum1_sum2, u_add_l_mul_u_add_l_sub_l_mul_l_u_mul_u, shift3);

    // release memory
    free(sum1_sum2.ptr_32);
//    if (shift3) free(sum3.type.ptr_64);
//    if (mul_upper_shift * 2) free(sum2.type.ptr_64);
    free(u_add_l_mul_u_add_l_sub_l_mul_l_u_mul_u.ptr_32);
    free(u_add_l_mul_u_add_l.ptr_32);
    free(subtrahend.ptr_32);
    free(u_add_l.ptr_32);
    free(l_mul_l.ptr_32);
    free(u_mul_u.ptr_32);

    return sum1_sum2_sum3;
}


/**
 * Multiply two 2*2 matrices. Each matrix has (n+1)-th, n-th,(n-1)-th elements.
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
        nxt_mul_nxt = karatubaMulNumber(a.f_n_nxt, b.f_n_nxt);
        n_mul_n = karatubaMulNumber(a.f_n_curr, b.f_n_curr);
        pre_mul_pre = karatubaMulNumber(a.f_n_prv, b.f_n_prv);
        m.f_n_nxt = addNumber(nxt_mul_nxt, n_mul_n, 0);
        m.f_n_prv = addNumber(n_mul_n, pre_mul_pre, 0);
        m.f_n_curr = subNumber(m.f_n_nxt, m.f_n_prv);
    } else {
        // square matrix multiplication
        nxt_mul_nxt = karatsubaMulSameNumber(a.f_n_nxt);
        n_mul_n = karatsubaMulSameNumber(a.f_n_curr);
        pre_mul_pre = karatsubaMulSameNumber(a.f_n_prv);
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
        printf("Couldn't allocate memory");
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
