#include "SWAVX_GPT_Opt.h"

#ifdef V1
// V1
void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal){
    int i, j;

    __m256i gapScoreVec = _mm256_set1_epi32(gapScore);
    __m256i matchScoreVec = _mm256_set1_epi32(matchScore);
    __m256i missmatchScoreVec = _mm256_set1_epi32(missmatchScore);

    for (i = 1; i < n; i++) {
        for (j = 1; j < m; j += 8) {
            __m256i upVec = _mm256_loadu_si256((__m256i*)&H[(i-1) * m + j]);
            upVec = _mm256_add_epi32(upVec, gapScoreVec);

            __m256i leftVec = _mm256_loadu_si256((__m256i*)&H[i * m + j - 1]);
            leftVec = _mm256_add_epi32(leftVec, gapScoreVec);

            __m256i diagVec = _mm256_loadu_si256((__m256i*)&H[(i-1) * m + j - 1]);
            __m256i aVec = _mm256_loadu_si256((__m256i*)&a[j-1]);
            __m256i bVec = _mm256_set1_epi8(b[i-1]);

            __m256i matchMissmatchVec = _mm256_cmpeq_epi8(aVec, bVec);
            matchMissmatchVec = _mm256_blendv_epi8(missmatchScoreVec, matchScoreVec, matchMissmatchVec);

            diagVec = _mm256_add_epi32(diagVec, matchMissmatchVec);

            __m256i maxVec = _mm256_max_epi32(upVec, leftVec);
            maxVec = _mm256_max_epi32(maxVec, diagVec);

            _mm256_storeu_si256((__m256i*)&H[i * m + j], maxVec);

            __m256i currentMaxVal = _mm256_loadu_si256((__m256i*)maxVal);
            __m256i newMaxVal = _mm256_max_epi32(currentMaxVal, maxVec);
            _mm256_storeu_si256((__m256i*)maxVal, newMaxVal);
        }
    }
}

int matchMissmatchScore(int i, int j, int8_t *a, int8_t *b) {
    return a[j-1] == b[i-1] ? matchScore : missmatchScore;
}

#endif