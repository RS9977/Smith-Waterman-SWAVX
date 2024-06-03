#include "SWAVX_Serial.h"

void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal){


    int i, j;
    
    for (i = 1; i < n; i++) { 
        for (j = 1; j < m; j++) { //
            similarityScore(a, b, m, n, i, j, H, P, maxVal);
        }
    }

}

void similarityScore(int8_t* a, int8_t* b, int m, int n, int i, int j, INT* H, INT* P, int* maxVal) {

    int up, left, diag;

    int index = m * i + j;

    up = H[index - m] + gapScore;

    left = H[index - 1] + gapScore;

    diag = H[index - m - 1] + matchMissmatchScore(i, j, a, b);


    int max = NONE;
    int pred = NONE;

    
    if (diag > max) { 
        max = diag;
        pred = DIAGONAL;
    }

    if (up > max) { 
        max = up;
        pred = UP;
    }
    
    if (left > max) {
        max = left;
        pred = LEFT;
    }
    H[index] = max;

    if (max > *maxVal) {
        *maxVal = H[index];
    }

}  

int matchMissmatchScore(int i, int j, int8_t*a, int8_t*b) {
    if (a[j-1] == b[i-1])
        return matchScore;
    else
        return missmatchScore;
} 