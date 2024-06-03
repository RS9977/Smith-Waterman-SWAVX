#include "SWAVX_GPT_Opt.h"

#ifdef V1
// V1 with the same name in functions / not correct
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

//V2 with the same name in functions / corrected version of V2
#elif V2
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal){
    int i, j;
    *maxVal = 0;

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
            __m256i aVec = _mm256_cvtepi8_epi32(_mm_loadl_epi64((__m128i*)&a[j-1]));
            __m256i bVec = _mm256_set1_epi32(b[i-1]);

            __m256i matchMissmatchVec = _mm256_cmpeq_epi32(aVec, bVec);
            matchMissmatchVec = _mm256_blendv_epi8(missmatchScoreVec, matchScoreVec, matchMissmatchVec);

            diagVec = _mm256_add_epi32(diagVec, matchMissmatchVec);

            __m256i maxVec = _mm256_max_epi32(upVec, leftVec);
            maxVec = _mm256_max_epi32(maxVec, diagVec);

            _mm256_storeu_si256((__m256i*)&H[i * m + j], maxVec);

            for (int k = 0; k < 8; k++) {
                int idx = i * m + j + k;
                if (H[idx] > *maxVal) {
                    *maxVal = H[idx];
                }
            }
        }
    }
}

int matchMissmatchScore(int i, int j, int8_t *a, int8_t *b) {
    return a[j-1] == b[i-1] ? matchScore : missmatchScore;
}


#elif V3

void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal){
    for (int i = 1; i < n; i++) { 
        for (int j = 1; j < m; j += 8) { 
            similarityScore(a, b, m, n, i, j, H, P, maxVal);
        }
    }
}

void similarityScore(int8_t* a, int8_t* b, int m, int n, int i, int j, INT* H, INT* P, int* maxVal) {
    int index = m * i + j;

    __m256i uu = _mm256_loadu_si256((__m256i*)&H[index - m]);
    __m256i ll = _mm256_loadu_si256((__m256i*)&H[index - 1]);
    __m256i dd = _mm256_set1_epi32(matchMissmatchScore(i, j, a, b));

    __m256i C_vec = _mm256_set1_epi32(gapScore);
    __m256i mm = _mm256_set1_epi32(NONE);

    uu = _mm256_add_epi32(uu, C_vec);
    ll = _mm256_add_epi32(ll, C_vec);
    dd = _mm256_add_epi32(dd, _mm256_loadu_si256((__m256i*)&H[index - m - 1]));

    __m256i pred = _mm256_set1_epi32(NONE);

    __m256i mask_dd = _mm256_cmpgt_epi32(dd, mm);
    mm = _mm256_blendv_epi8(mm, dd, mask_dd);
    pred = _mm256_blendv_epi8(pred, _mm256_set1_epi32(DIAGONAL), mask_dd);

    __m256i mask_uu = _mm256_cmpgt_epi32(uu, mm);
    mm = _mm256_blendv_epi8(mm, uu, mask_uu);
    pred = _mm256_blendv_epi8(pred, _mm256_set1_epi32(UP), mask_uu);

    __m256i mask_ll = _mm256_cmpgt_epi32(ll, mm);
    mm = _mm256_blendv_epi8(mm, ll, mask_ll);
    pred = _mm256_blendv_epi8(pred, _mm256_set1_epi32(LEFT), mask_ll);

    _mm256_storeu_si256((__m256i*)&H[index], mm);

    __m256i maxVal_vec = _mm256_loadu_si256((__m256i*)maxVal);
    __m256i mask_max = _mm256_cmpgt_epi32(mm, maxVal_vec);
    maxVal_vec = _mm256_blendv_epi8(maxVal_vec, mm, mask_max);
    _mm256_storeu_si256((__m256i*)maxVal, maxVal_vec);
}

int matchMissmatchScore(int i, int j, int8_t*a, int8_t*b) {
    return (a[j-1] == b[i-1]) ? matchScore : missmatchScore;
}

//V4 asking for corrected version of the V3

#elif V4

void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal) {
    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m; j += 8) {  // Process 8 elements at a time
            similarityScore(a, b, m, n, i, j, H, P, maxVal);
        }
    }
}

void similarityScore(int8_t *a, int8_t *b, int m, int n, int i, int j, INT *H, INT *P, int *maxVal) {
    int index = m * i + j;

    // Load 8 elements into AVX2 registers
    __m256i uu = _mm256_loadu_si256((__m256i*)&H[index - m]);
    __m256i ll = _mm256_loadu_si256((__m256i*)&H[index - 1]);
    __m256i dd = _mm256_set1_epi32(matchMissmatchScore(i, j, a, b));

    __m256i C_vec = _mm256_set1_epi32(gapScore);
    __m256i mm = _mm256_set1_epi32(NONE);

    uu = _mm256_add_epi32(uu, C_vec);
    ll = _mm256_add_epi32(ll, C_vec);
    dd = _mm256_add_epi32(dd, _mm256_loadu_si256((__m256i*)&H[index - m - 1]));

    __m256i pred = _mm256_set1_epi32(NONE);

    __m256i mask_dd = _mm256_cmpgt_epi32(dd, mm);
    mm = _mm256_blendv_epi8(mm, dd, mask_dd);
    pred = _mm256_blendv_epi8(pred, _mm256_set1_epi32(DIAGONAL), mask_dd);

    __m256i mask_uu = _mm256_cmpgt_epi32(uu, mm);
    mm = _mm256_blendv_epi8(mm, uu, mask_uu);
    pred = _mm256_blendv_epi8(pred, _mm256_set1_epi32(UP), mask_uu);

    __m256i mask_ll = _mm256_cmpgt_epi32(ll, mm);
    mm = _mm256_blendv_epi8(mm, ll, mask_ll);
    pred = _mm256_blendv_epi8(pred, _mm256_set1_epi32(LEFT), mask_ll);

    _mm256_storeu_si256((__m256i*)&H[index], mm);

    __m256i maxVal_vec = _mm256_loadu_si256((__m256i*)maxVal);
    __m256i mask_max = _mm256_cmpgt_epi32(mm, maxVal_vec);
    maxVal_vec = _mm256_blendv_epi8(maxVal_vec, mm, mask_max);
    _mm256_storeu_si256((__m256i*)maxVal, maxVal_vec);
}

int matchMissmatchScore(int i, int j, int8_t *a, int8_t *b) {
    return (a[j - 1] == b[i - 1]) ? matchScore : missmatchScore;
}


//use Compiler Sage of Chat GPT
#elif V5
void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal) {
    int i, j;

    for (i = 1; i < n; i++) {
        for (j = 1; j < m; j+=8) {
            similarityScore(a, b, m, n, i, j, H, P, maxVal);
        }
    }
}

void similarityScore(int8_t* a, int8_t* b, int m, int n, int i, int j, INT* H, INT* P, int* maxVal) {
    __m256i C_vec = _mm256_set1_epi32(gapScore);
    __m256i max_val_vec = _mm256_set1_epi32(*maxVal);
    __m256i none_vec = _mm256_set1_epi32(NONE);

    for (int k = 0; k < 8; k++) {
        int index = m * i + j + k;

        __m256i H_vec = _mm256_loadu_si256((__m256i*)&H[index - m]);
        __m256i H_left_vec = _mm256_loadu_si256((__m256i*)&H[index - 1]);

        __m256i uu_vec = _mm256_add_epi32(H_vec, C_vec);
        __m256i ll_vec = _mm256_add_epi32(H_left_vec, C_vec);

        int function3_val = matchMissmatchScore(i, j + k, a, b);
        __m256i dd_vec = _mm256_add_epi32(H_vec, _mm256_set1_epi32(function3_val));

        __m256i mm_vec = none_vec;
        __m256i pred_vec = none_vec;

        __m256i mask_dd = _mm256_cmpgt_epi32(dd_vec, mm_vec);
        mm_vec = _mm256_blendv_epi8(mm_vec, dd_vec, mask_dd);
        pred_vec = _mm256_blendv_epi8(pred_vec, _mm256_set1_epi32(DIAGONAL), mask_dd);

        __m256i mask_uu = _mm256_cmpgt_epi32(uu_vec, mm_vec);
        mm_vec = _mm256_blendv_epi8(mm_vec, uu_vec, mask_uu);
        pred_vec = _mm256_blendv_epi8(pred_vec, _mm256_set1_epi32(UP), mask_uu);

        __m256i mask_ll = _mm256_cmpgt_epi32(ll_vec, mm_vec);
        mm_vec = _mm256_blendv_epi8(mm_vec, ll_vec, mask_ll);
        pred_vec = _mm256_blendv_epi8(pred_vec, _mm256_set1_epi32(LEFT), mask_ll);

        _mm256_storeu_si256((__m256i*)&H[index], mm_vec);

        __m256i mask_max = _mm256_cmpgt_epi32(mm_vec, max_val_vec);
        max_val_vec = _mm256_blendv_epi8(max_val_vec, mm_vec, mask_max);
    }

    __m256i temp = _mm256_permute2x128_si256(max_val_vec, max_val_vec, 1);
    max_val_vec = _mm256_max_epi32(max_val_vec, temp);
    temp = _mm256_shuffle_epi32(max_val_vec, _MM_SHUFFLE(2, 3, 0, 1));
    max_val_vec = _mm256_max_epi32(max_val_vec, temp);
    temp = _mm256_shuffle_epi32(max_val_vec, _MM_SHUFFLE(1, 0, 3, 2));
    max_val_vec = _mm256_max_epi32(max_val_vec, temp);

    *maxVal = _mm256_extract_epi32(max_val_vec, 0);
}

int matchMissmatchScore(int i, int j, int8_t* a, int8_t* b) {
    return (a[j-1] == b[i-1]) ? matchScore : missmatchScore;
}

#endif