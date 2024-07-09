cmake_minimum_required(VERSION 3.25)
project(CalcFibonacci C)

set(CMAKE_C_STANDARD 17)

add_executable(CalcFibonacci main.c calcFibonacci.c calcFibonacci.h util.c util.h bigInt.h)
