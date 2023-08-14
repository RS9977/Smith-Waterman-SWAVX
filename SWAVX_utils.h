#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

#define RESET   "\033[0m"
#define BOLDRED "\033[1m\033[31m"      /* Bold Red */

void printMatrix(int* matrix, int8_t *a, int8_t *b, int m, int n);
void printPredecessorMatrix(int* matrix, int8_t *a, int8_t *b, int m, int n);
void generate(int8_t *a, int8_t *b, int m, int n);