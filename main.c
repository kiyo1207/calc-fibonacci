#include <stdbool.h>
#include <stddef.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include "calcFibonacci.h"
#include "util.h"


const char *usage_msg =
        "Usage: %s [options] -n X   Calculate the n-th Fibonacci number.\n"
        "   or: %s -h               Show help message and exit.\n"
        "   or: %s --help           Same as '-h' option.\n";

const char *help_msg =
        "Optional arguments:\n"
        "  -n X    Position of the fibonacci number to be calculated.\n"
        "          The first fibonacci number (n=0) is 0, the second (n=1) is 1. The maximum number of n is 2^64-1.\n"
        "          The Number X can be given in decimal, octal, and hexadecimal.\n"
        "  -d X    Display calculation results(default: X = 0).\n"
        "          0:  Do not display.\n"
        "          1:  In hexadecimal.\n"
        "          2:  In decimal(The larger the number, the longer it takes).\n"
        "  -t      Measure the runtime of the calculation. The measured time is from the beginning of the calculation to"
        "          writing the result in memory. The tun time for displaying the calculation result is not included."
        "  -h      Show help message (this text) and exit.\n"
        "  --help  Same as '-h' option.\n";

void print_usage(const char *progname) {
    fprintf(stderr, usage_msg, progname, progname, progname);
}

void print_help(const char *progname) {
    print_usage(progname);
    fprintf(stderr, "\n%s", help_msg);
}


int main(int argc, char **argv) {

    const char *programmName = argv[0];
    if (argc == 1) {
        print_usage(programmName);
        return EXIT_FAILURE;
    }

    struct option longOpts[] = {
            {"help", no_argument, NULL, 'h'},
            {NULL, 0,             NULL, 0}
    };

    uint64_t n;
    int display = 0;
    bool benchmarkFlg = false;
    bool nSetFlg = false;
    int opt;
    int longIndex;


    while ((opt = getopt_long(argc, argv, "n:d:th", longOpts, &longIndex)) != -1) {

        switch (opt) {
            case 'n':
                errno = 0;
                char *endPtrN;
                n = strtoull(optarg, &endPtrN, 0);
                if (*endPtrN != 0) {
                    fprintf(stderr, "Invalid number for '-n' -- %s\n", optarg);
                    print_usage(programmName);
                    return EXIT_FAILURE;
                }

                if (errno != 0 || optarg[0] == '-') {
                    fprintf(stderr,
                            "The number of fibonacci must be between 0 and 2^64-1.\n");
                    print_usage(programmName);
                    return EXIT_FAILURE;
                }
                nSetFlg = true;
                break;

            case 'd':
                if (strlen(optarg) == 1 && isdigit(optarg[0])) {
                    display = atoi(optarg);
                    if (display == 0 || display == 1 || display == 2) {
                        break;
                    }
                }
                fprintf(stderr, "Invalid result's display mode for '-d' -- %s\n", optarg);
                print_usage(programmName);
                return EXIT_FAILURE;

            case 't':
                benchmarkFlg = true;
                break;

            case 'h':
                print_help(programmName);
                return EXIT_SUCCESS;

            default : /* ’? ’ */
                print_usage(programmName);
                return EXIT_FAILURE;
        }
    }


    if (!nSetFlg) {
        fprintf(stderr, "The number of fibonacci for '-n' is missing.\n");
        print_usage(programmName);
        return EXIT_FAILURE;
    }


    // Approximation of 256^(k-1) < (1/sqrt(5) * ((1 + sqrt(5)) / 2)^n for k.
    size_t size = n / 11 + 1;
    uint8_t *buffer = malloc(size * sizeof(uint8_t));
    if (buffer == NULL) {
        printf("Couldn't allocate memory : main().\n");
        abort();
    }
    struct timespec start, end;
    double time;
    if (benchmarkFlg) {
        clock_gettime(CLOCK_MONOTONIC, &start);
        fib(n, size, buffer);
        clock_gettime(CLOCK_MONOTONIC, &end);
        time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);
        printf("Calculation finished.[time: %lfs]\n", time);
    } else {
        fib(n, size, buffer);
        printf("Calculation finished.\n");
    }


    if (display == 1) {
        printf("Result:0x%s\n", uint8BufferToHexString(size, buffer));
    } else if (display == 2) {
        printf("Result:%s\n", uint8BufferToDecString(size, buffer));
    }

    return EXIT_SUCCESS;

}