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
        printf("Given char is no hex digit: %c\n", letter);
        abort();
    }
}


char *uint8BufferToHexString(size_t n, const uint8_t *ptr) {
    //Reserve memory for the return value
    char *ret = malloc(n * 2 * sizeof(char) + 1);
    if (ret == NULL) {
        printf("Could not reserve memory for HexString");
        abort();
    }
    bool firstdigitWritten = false; //shortens leading zeros
    size_t stringIndex = 0;

    //starts at the most significant bit and converts one u_int8t value to two hex chars
    for (size_t i = 0; i < n; i++) {
        uint8_t val = *(ptr + n - 1 - i);

        if ((val >> 4) > 0 || firstdigitWritten) {
            *(ret + stringIndex++) = intToHex(val >> 4);
            firstdigitWritten = true;
        }
        if ((val & 0xf) > 0 || firstdigitWritten) {
            *(ret + stringIndex++) = intToHex(val & 0xf);
            firstdigitWritten = true;
        }
    }

    // return zero if no digits were written
    if (!firstdigitWritten) {
        *(ret + stringIndex++) = '0';
    }

    //set the last arrayfield to zero for null termination of String
    *(ret + stringIndex) = 0;
    return ret;
}

char *bigIntToHexString(struct bigInt *a) {
    return uint8BufferToHexString(a->size_32 * 4, (uint8_t *) a->ptr_32);
}



struct bigInt hexStringToBigInt(const char *string) {
    size_t stringLength = strlen(string);
    //leerer string
    if (stringLength == 0) {
        printf("Tried to create a bigInt from emptyString\n");
        abort();
    }

    size_t arraylength = (stringLength) / 8;
    if (stringLength % 8 > 0) arraylength++;


    //reserve memory
    uint8_t *array = calloc(arraylength, sizeof(uint32_t));
    if (array == 0) {
        printf("Could'nt allocate memory for bigInt.\n");
        abort();
    }

    size_t stringIndex = stringLength;
    size_t arrayIndex = 0;
    uint8_t buffer = 0;
    bool bufferused = false;

    //goes through the string from tail (=least significant bits) and converts two chars to one u_int8_t number
    do {
        stringIndex--;
        char letter = *(string + stringIndex);

        if (bufferused) {
            buffer += hexToInt(letter) << 4;
            bufferused = false;
            *(array + arrayIndex) = buffer;
            arrayIndex++;
        } else {
            buffer = hexToInt(letter);
            bufferused = true;
        }
    } while (stringIndex != 0);

    if (bufferused) *(array + arrayIndex) = buffer;

    return (struct bigInt) {.size_32 = arraylength, .ptr_32 = (uint32_t *) array};
}

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
                    fprintf(stderr, "Could'nt allocate memory for the decimal string");
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

    struct bigInt quotient;
    quotient.size_32 = (n / 4) + (n % 4 == 0 ? 0 : 1);
    quotient.ptr_32 = calloc(b.size_32, sizeof(uint32_t));

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

    // if 0th elemente(lowest position) and the result of division is 0, terminate the process
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