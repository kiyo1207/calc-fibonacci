#ifndef CALCFIBONACCI_UTIL_H
#define CALCFIBONACCI_UTIL_H



#include "bigInt.h"

char intToHex(uint8_t figure);
uint8_t hexToInt(char letter);char *uint8BufferToHexString(size_t n, const uint8_t *ptr);
char *bigIntToHexString(struct bigInt *a);
struct bigInt hexStringToBigInt(const char *string);
char *uint8BufferToDecString(size_t n, const uint8_t *buf);

#endif //CALCFIBONACCI_UTIL_H
