#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include <stdint.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include <stdbool.h>
#include <sys/stat.h>

#include "SWAVX_utils.h"
#include "SWAVX_SubMat.h"

#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

void SWAVX_512_SeqToSeq_SubMat(int8_t *a, int8_t *b, int *H, int* P, int m, int n, int NumOfTest, int gapScore);
void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, int* H, int* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, int gapScore, int8_t *a, int8_t *b, int m, int n);
void similarityScoreIntrinsic(__m512i* HH,__m512i* Hu,__m512i* Hd,__m512i* Hl,__m512i* PP, __m512i reverseIndices, long long int ii, long long int jj, int* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, int gapScore, int8_t *a, int8_t *b, int m, int n);
void backtrack(int* P, long long int maxPos, long long int maxPos_max_len, int m, int n);
int matchMissmatchScore(long long int i, long long int j, int8_t* a, int8_t* b);