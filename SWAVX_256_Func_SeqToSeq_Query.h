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

#define BLOCKSIZE  4000

#ifdef SUBMAT
typedef __m128i RevType;
#else
typedef __m256i RevType;
#endif


_Thread_local extern __m256i scores_mat[32];
_Thread_local extern __m256i local_max;

void SWAVX_256_SeqToSeq_QueryBTBatch(int8_t *a, ProteinBatch batch, INT *H, INT* P, int m, int NumOfTest, INT* maxVal, int8_t* query_prof, int MaxHSize);
void SWAVX_256_SeqToSeq_QueryLTBatch(ProteinBatch batch, int8_t *b, INT *H, INT* P, int n, int NumOfTest, INT* maxVal, int8_t* query_prof, int MaxHSize);
void similarityScore_QueryBTBatch(int ind, int ind_u, int ind_d, int ind_l, int ii, int jj, INT* H, INT* P, int max_len, int* maxPos, int *maxPos_max_len, INT* maxVal, int8_t *a, ProteinBatch batch, int m, int k);
void similarityScore_QueryLTBatch(int ind, int ind_u, int ind_d, int ind_l, int ii, int jj, INT* H, INT* P, int max_len, int* maxPos, int *maxPos_max_len, INT* maxVal, ProteinBatch batch, int8_t* b, int n, int k);
void similarityScoreIntrinsic8(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, int ii, int jj, INT* H, int ind, int max_len, int* maxPos, int *maxPos_max_len,  INT* maxVal, int k, int score_bias, int reminder);
void backtrack(INT* P, int maxPos, int maxPos_max_len, int m, int n);
int matchMissmatchScore(int i, int j, int8_t* a, int8_t* b);
void Transpose32By32(__m256i mat[32]);