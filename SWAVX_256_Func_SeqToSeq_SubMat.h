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


#ifdef SUBMAT
typedef __m128i RevType;
#else
typedef __m256i RevType;
#endif

_Thread_local extern __m256i local_max;

void SWAVX_256_SeqToSeq_SubMat(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal, int8_t* query_prof, bool query_case);
void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, INT* H, INT* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, int m, int n);
void similarityScoreIntrinsic32(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, long long int ii, long long int jj, INT* H, long long int ind,  long long int max_len, long long int* maxPos,  long long int *maxPos_max_len,  INT* maxVal, int8_t *a, int8_t *b, int m, int n);
void similarityScoreIntrinsic16(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m128i reverseIndices, long long int ii, long long int jj, INT* H, long long int ind,  long long int max_len, long long int* maxPos,  long long int *maxPos_max_len,  INT* maxVal, int8_t *a, int8_t *b, int m, int n);
void similarityScoreIntrinsic8(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP,  int8_t* scores, RevType reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len,  INT* maxVal, int8_t *a, int8_t *b, int m, int n, int k, bool query_case);
void similarityScoreIntrinsic32_affine(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind,  long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n);
void similarityScoreIntrinsic16_affine(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m128i reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n, int k);
void similarityScoreIntrinsic8_affine(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, RevType reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n, int k);
void similarityScore_affine(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, INT* H, INT* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n);
void backtrack(INT* P, long long int maxPos, long long int maxPos_max_len, int m, int n);
int matchMissmatchScore(long long int i, long long int j, int8_t* a, int8_t* b);