#include "SWAVX_Serial.h"

/*--------------------------------------------------------------------
 * Function:    main
 */
void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal){



    //Calculates the similarity matrix
    int i, j;
    
    for (i = 1; i < n; i++) { //Lines
        for (j = 1; j < m; j++) { //Columns
            similarityScore(a, b, m, n, i, j, H, P, maxVal);
        }
    }

}


/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate  the maximum Similarity-Score H(i,j)
 */
void similarityScore(int8_t* a, int8_t* b, int m, int n, int i, int j, INT* H, INT* P, int* maxVal) {

    int up, left, diag;

    //Stores index of element
    int index = m * i + j;

    //Get element above
    up = H[index - m] + gapScore;

    //Get element on the left
    left = H[index - 1] + gapScore;

    //Get element on the diagonal
    diag = H[index - m - 1] + matchMissmatchScore(i, j, a, b);

    //Calculates the maximum
    int max = NONE;
    int pred = NONE;

    
    if (diag > max) { //same letter ↖
        max = diag;
        pred = DIAGONAL;
    }

    if (up > max) { //remove letter ↑ 
        max = up;
        pred = UP;
    }
    
    if (left > max) { //insert letter ←
        max = left;
        pred = LEFT;
    }
    //Inserts the value in the similarity and predecessor matrixes
    H[index] = max;

    //Updates maximum score to be used as seed on backtrack 
    if (max > *maxVal) {
        *maxVal = H[index];
    }

}  /* End of similarityScore */


/*--------------------------------------------------------------------
 * Function:    matchMissmatchScore
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
int matchMissmatchScore(int i, int j, int8_t*a, int8_t*b) {
    if (a[j-1] == b[i-1])
        return matchScore;
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */