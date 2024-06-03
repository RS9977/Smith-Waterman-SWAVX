#include "SWAVX_GPT.h"

#ifdef V1
// V1
void SWAVX_GPT_Func(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int *maxVal) {
    __m256i zero = _mm256_setzero_si256();
    __m256i maxScore = zero;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; j += 8) {
            __m256i score_diag = _mm256_loadu_si256((__m256i*)&H[(i - 1) * (n + 1) + (j - 1)]);
            __m256i score_up = _mm256_loadu_si256((__m256i*)&H[(i - 1) * (n + 1) + j]);
            __m256i score_left = _mm256_loadu_si256((__m256i*)&H[i * (n + 1) + (j - 1)]);

            __m256i a_elem = _mm256_set1_epi8(a[i - 1]);
            __m256i b_elems = _mm256_loadu_si256((__m256i*)&b[j - 1]);
            __m256i match_mismatch = _mm256_cmpeq_epi8(a_elem, b_elems);
            __m256i match = _mm256_blendv_epi8(_mm256_set1_epi8(missmatchScore), _mm256_set1_epi8(matchScore), match_mismatch);

            score_diag = _mm256_add_epi32(score_diag, _mm256_cvtepi8_epi32(_mm256_extracti128_si256(match, 0)));
            score_up = _mm256_add_epi32(score_up, _mm256_set1_epi32(gapScore));
            score_left = _mm256_add_epi32(score_left, _mm256_set1_epi32(gapScore));

            __m256i result = _mm256_max_epi32(zero, _mm256_max_epi32(score_diag, _mm256_max_epi32(score_up, score_left)));
            _mm256_storeu_si256((__m256i*)&H[i * (n + 1) + j], result);

            maxScore = _mm256_max_epi32(maxScore, result);
        }
    }

    // Find the maximum score in the maxScore vector
    int maxScoreArray[8];
    _mm256_storeu_si256((__m256i*)maxScoreArray, maxScore);
    for (int i = 0; i < 8; ++i) {
        if (maxScoreArray[i] > *maxVal) {
            *maxVal = maxScoreArray[i];
        }
    }
}


#elif V2
//V2 Optimizaed, incorrect
void SWAVX_GPT_Func(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int *maxVal) {
     __m256i zero = _mm256_setzero_si256();
    __m256i maxScore = zero;
    
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; j += 8) {
            // Load the values from the scoring matrix
            __m256i score_diag = _mm256_loadu_si256((__m256i*)&H[(i - 1) * (n + 1) + (j - 1)]);
            __m256i score_up = _mm256_loadu_si256((__m256i*)&H[(i - 1) * (n + 1) + j]);
            __m256i score_left = _mm256_loadu_si256((__m256i*)&H[i * (n + 1) + (j - 1)]);

            // Load the current characters of the sequences
            __m256i a_elem = _mm256_set1_epi8(a[i - 1]);
            __m256i b_elems = _mm256_loadu_si256((__m256i*)&b[j - 1]);

            // Calculate match/mismatch scores
            __m256i match_mismatch = _mm256_cmpeq_epi8(a_elem, b_elems);
            __m256i match = _mm256_blendv_epi8(_mm256_set1_epi8(missmatchScore), _mm256_set1_epi8(matchScore), match_mismatch);

            // Convert the int8_t match results to int32_t
            __m256i match_low = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(match, 0));
            __m256i match_high = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(match, 1));

            // Calculate scores
            __m256i result_low = _mm256_add_epi32(score_diag, match_low);
            __m256i result_high = _mm256_add_epi32(score_diag, match_high);
            
            __m256i gap_score = _mm256_set1_epi32(gapScore);
            result_low = _mm256_max_epi32(result_low, _mm256_add_epi32(score_up, gap_score));
            result_low = _mm256_max_epi32(result_low, _mm256_add_epi32(score_left, gap_score));
            result_low = _mm256_max_epi32(result_low, zero);
            
            result_high = _mm256_max_epi32(result_high, _mm256_add_epi32(score_up, gap_score));
            result_high = _mm256_max_epi32(result_high, _mm256_add_epi32(score_left, gap_score));
            result_high = _mm256_max_epi32(result_high, zero);

            // Store the results back to the scoring matrix
            _mm256_storeu_si256((__m256i*)&H[i * (n + 1) + j], result_low);
            _mm256_storeu_si256((__m256i*)&H[i * (n + 1) + j + 4], result_high);

            // Update the maximum score
            maxScore = _mm256_max_epi32(maxScore, result_low);
            maxScore = _mm256_max_epi32(maxScore, result_high);
        }
    }

    // Extract the maximum score from the maxScore vector
    int maxScoreArray[8];
    _mm256_storeu_si256((__m256i*)maxScoreArray, maxScore);
    for (int i = 0; i < 8; ++i) {
        if (maxScoreArray[i] > *maxVal) {
            *maxVal = maxScoreArray[i];
        }
    }
}

#elif V3
//V3 optimized
void SWAVX_GPT_Func(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int *maxVal) {
    __m256i zero = _mm256_setzero_si256();
    __m256i maxScore = zero;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; j += 8) {
            // Load the values from the scoring matrix
            __m256i score_diag = _mm256_loadu_si256((__m256i*)&H[(i - 1) * (n + 1) + (j - 1)]);
            __m256i score_up = _mm256_loadu_si256((__m256i*)&H[(i - 1) * (n + 1) + j]);
            __m256i score_left = _mm256_loadu_si256((__m256i*)&H[i * (n + 1) + (j - 1)]);

            // Load the current characters of the sequences
            __m256i a_elem = _mm256_set1_epi8(a[i - 1]);
            __m256i b_elems = _mm256_loadu_si256((__m256i*)&b[j - 1]);

            // Calculate match/mismatch scores
            __m256i match_mismatch = _mm256_cmpeq_epi8(a_elem, b_elems);
            __m256i match = _mm256_blendv_epi8(_mm256_set1_epi8(missmatchScore), _mm256_set1_epi8(matchScore), match_mismatch);

            // Convert the int8_t match results to int32_t
            __m256i match_low = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(match, 0));
            __m256i match_high = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(match, 1));

            // Split scores for lower and upper halves of the vectors
            __m256i score_diag_low = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(score_diag, 0));
            __m256i score_diag_high = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(score_diag, 1));
            __m256i score_up_low = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(score_up, 0));
            __m256i score_up_high = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(score_up, 1));
            __m256i score_left_low = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(score_left, 0));
            __m256i score_left_high = _mm256_cvtepi8_epi32(_mm256_extracti128_si256(score_left, 1));

            // Calculate scores
            __m256i result_low = _mm256_add_epi32(score_diag_low, match_low);
            __m256i result_high = _mm256_add_epi32(score_diag_high, match_high);
            
            __m256i gap_score = _mm256_set1_epi32(gapScore);
            result_low = _mm256_max_epi32(result_low, _mm256_add_epi32(score_up_low, gap_score));
            result_low = _mm256_max_epi32(result_low, _mm256_add_epi32(score_left_low, gap_score));
            result_low = _mm256_max_epi32(result_low, zero);
            
            result_high = _mm256_max_epi32(result_high, _mm256_add_epi32(score_up_high, gap_score));
            result_high = _mm256_max_epi32(result_high, _mm256_add_epi32(score_left_high, gap_score));
            result_high = _mm256_max_epi32(result_high, zero);

            // Store the results back to the scoring matrix
            _mm256_storeu_si256((__m256i*)&H[i * (n + 1) + j], result_low);
            _mm256_storeu_si256((__m256i*)&H[i * (n + 1) + j + 4], result_high);

            // Update the maximum score
            maxScore = _mm256_max_epi32(maxScore, result_low);
            maxScore = _mm256_max_epi32(maxScore, result_high);
        }
    }

    // Extract the maximum score from the maxScore vector
    int maxScoreArray[8];
    _mm256_storeu_si256((__m256i*)maxScoreArray, maxScore);
    for (int i = 0; i < 8; ++i) {
        if (maxScoreArray[i] > *maxVal) {
            *maxVal = maxScoreArray[i];
        }
    }

}

// One prompt
#elif V4
void SWAVX_GPT_Func(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int *maxVal) {
     __m256i vGap = _mm256_set1_epi32(gapScore);
    __m256i vZero = _mm256_setzero_si256();
    __m256i vMaxVal = vZero;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; j += 8) {
            // Load sequences and previous scores
            __m256i vH_left = _mm256_loadu_si256((__m256i*)&H[(i-1) * (n+1) + j]);
            __m256i vH_up = _mm256_loadu_si256((__m256i*)&H[i * (n+1) + (j-1)]);

            // Load characters
            __m256i vA = _mm256_set1_epi32(a[i-1]);
            __m256i vB = _mm256_loadu_si256((__m256i*)&b[j-1]);

            // Compute match/mismatch scores
            __m256i vMatchMask = _mm256_cmpeq_epi8(vA, vB);
            __m256i vMismatchMask = _mm256_xor_si256(vMatchMask, _mm256_set1_epi8(-1));

            __m256i vMatchScore = _mm256_set1_epi32(matchScore);
            __m256i vMismatchScore = _mm256_set1_epi32(missmatchScore);

            __m256i vScores = _mm256_blendv_epi8(vMismatchScore, vMatchScore, vMatchMask);

            // Update scores
            vH_left = _mm256_add_epi32(vH_left, vGap);
            vH_up = _mm256_add_epi32(vH_up, vGap);

            __m256i vDiag = _mm256_loadu_si256((__m256i*)&H[(i-1) * (n+1) + (j-1)]);
            vDiag = _mm256_add_epi32(vDiag, vScores);

            __m256i vScore = _mm256_max_epi32(vDiag, vZero);
            vScore = _mm256_max_epi32(vScore, vH_left);
            vScore = _mm256_max_epi32(vScore, vH_up);

            // Store scores
            _mm256_storeu_si256((__m256i*)&H[i * (n+1) + j], vScore);

            // Update max value
            vMaxVal = _mm256_max_epi32(vMaxVal, vScore);
        }
    }

    // Extract max value from vMaxVal
    for (int k = 0; k < 8; ++k) {
        int val = _mm256_extract_epi32(vMaxVal, k);
        if (val > *maxVal) {
            *maxVal = val;
        }
    }
}


#endif