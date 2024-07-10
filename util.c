#include <stdint.h>
#include <printf.h>
#include <stdlib.h>
#include <stdbool.h>
#include <libc.h>
#include <math.h>
#include "util.h"


char intToHex(uint8_t figure) {
    if (figure < 10) {
        return figure + 48;
    } else {
        return figure - 10 + 65;
    }
}


uint8_t hexToInt(char letter) {
    if ('0' <= letter && letter <= '9') {
        return letter - 48;
    } else if ('a' <= letter && letter <= 'f') {
        return letter - 97 + 10;
    } else if ('A' <= letter && letter <= 'F') {
        return letter - 65 + 10;
    } else {
        printf("Given char is not valid hex character.-- %c\n", letter);
        abort();
    }
}


/**
 *
 * @param n
 * @param ptr
 * @return
 */
char *uint8BufferToHexString(size_t n, const uint8_t *ptr) {

    char *ret = malloc(n * 2 * sizeof(char) + 1);
    if (ret == NULL) {
        printf("Couldn't allocate memory : uint8BufferToHexString().\n");
        abort();
    }
    bool notZeroFlg = false; //shortens leading zeros
    size_t stringIndex = 0;

    //starts at the most significant bit and converts one u_int8t value to two hex chars
    for (size_t i = 0; i < n; i++) {
        uint8_t val = *(ptr + n - 1 - i);

        if ((val >> 4) > 0 || notZeroFlg) {
            *(ret + stringIndex++) = intToHex(val >> 4);
            notZeroFlg = true;
        }
        if ((val & 0xf) > 0 || notZeroFlg) {
            *(ret + stringIndex++) = intToHex(val & 0xf);
            notZeroFlg = true;
        }
    }

    // return zero if no digits were written
    if (!notZeroFlg) {
        *(ret + stringIndex++) = '0';
    }

    //set the last arrayfield to zero for null termination of String
    *(ret + stringIndex) = 0;
    return ret;

}

/**
 *
 * @param n : Length of buffer.
 * @param buf : Memory in which big integer is stored in little-endian.
 * @return big integer in decimal notation.
 */
char *uint8BufferToDecString(size_t n, const uint8_t *buf) {

    uint32_t divider = pow(10.0, 9);
    uint64_t remainder = 0;

    // Leading Zero padding is removed.
    uint64_t zeroPaddingCnt = 0;
    for (size_t i = n - 1; true; i--) {
        if (buf[i] == 0) {
            zeroPaddingCnt++;
            if (i != 0) {
                continue;
            } else {
                char *ch = calloc(1, sizeof(char));
                if (ch == NULL) {
                    printf("Couldn't allocate memory : uint8BufferToDecString().\n");
                    abort();
                }
                ch[0] = '0';
                return ch;
            }
        }
        break;
    }
    n -= zeroPaddingCnt;

    struct bigInt b;
    b.size_32 = (n / 4) + (n % 4 == 0 ? 0 : 1);
    b.ptr_32 = calloc(b.size_32, sizeof(uint32_t));
    if (b.ptr_32 == NULL) {
        printf("Couldn't allocate memory : uint8BufferToDecString().\n");
        abort();
    }

    struct bigInt quotient;
    quotient.size_32 = (n / 4) + (n % 4 == 0 ? 0 : 1);
    quotient.ptr_32 = calloc(b.size_32, sizeof(uint32_t));
    if (quotient.ptr_32 == NULL) {
        printf("Couldn't allocate memory : uint8BufferToDecString().\n");
        abort();
    }

    size_t k = 0;
    for (; k + 3 < n;) {
        b.ptr_32[k / 4] = ((uint32_t) buf[k]) + (((uint32_t) buf[k + 1]) << 8) + (((uint32_t) buf[k + 2]) << 16) +
                          (((uint32_t) buf[k + 3]) << 24);
        k += 4;
    }

    int rest = n % 4;
    for (int j = rest - 1; j >= 0; j--) {
        b.ptr_32[b.size_32 - 1] = (b.ptr_32[b.size_32 - 1] << 8) + (uint32_t) buf[k + j];
    }

    // array of decimal integer
    size_t digitDecimal = (unsigned long) ceil(n * 8 * log10(2));
    char *decimalCharArray = calloc(digitDecimal, sizeof(char));
    size_t decimalArrayCnt = 0;

    // if 0th element(the least significant position) and the result of division is 0, terminate the process
    while (!(b.ptr_32[b.size_32 - 1] / divider == 0 && b.size_32 - 1 == 0)) {

        int lag = 0;
        if (((b.ptr_32[b.size_32 - 1]) / divider) == 0) {
            quotient.size_32--;
            remainder = b.ptr_32[b.size_32 - 1];
            lag++;
        }
        for (int j = 0; j < b.size_32 - lag; j++) {
            quotient.ptr_32[quotient.size_32 - (j + 1)] =
                    ((remainder << 32) + b.ptr_32[b.size_32 - (j + lag + 1)]) / divider;
            remainder = ((uint64_t) (remainder << 32) + (uint64_t) b.ptr_32[b.size_32 - (j + lag + 1)]) -
                        ((uint64_t) quotient.ptr_32[quotient.size_32 - (j + 1)] * (uint64_t) divider);
        }


        // quotient becomes the next next dividend.
        b.size_32 = quotient.size_32;
        for (int i = 0; i < b.size_32; i++) {
            b.ptr_32[i] = quotient.ptr_32[i];
            quotient.ptr_32[i] = 0;
        }

        //　convert 10^9 decimal digits to characters.
        for (int i = 0; i < 9; i++) {
            decimalCharArray[decimalArrayCnt + i] = (remainder % 10) + '0';
            remainder /= 10;
        }
        decimalArrayCnt += 9;
        remainder = 0;
    }

    if (b.ptr_32[b.size_32 - 1] != 0) {
        uint32_t tmp = b.ptr_32[b.size_32 - 1];
        while (tmp) {
            decimalCharArray[decimalArrayCnt++] = (tmp % 10) + '0';
            tmp /= 10;
        }
    }
    for (size_t i = 0, j = decimalArrayCnt - 1; i < j; i++, j--) {
        char tmp = decimalCharArray[i];
        decimalCharArray[i] = decimalCharArray[j];
        decimalCharArray[j] = tmp;
    }

    free(quotient.ptr_32);
    free(b.ptr_32);
    return decimalCharArray;
}