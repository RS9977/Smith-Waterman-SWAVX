#include "SWAVX_256_Func_SeqToSeq_SubMat.h"

void SWAVX_256_SeqToSeq_SubMat(int8_t *a, int8_t *b, int *H, int* P, int m, int n, int NumOfTest, int gapScore){

    //Calculates the similarity matrix
    long long int i, j;

    //Start position for backtrack
    #ifdef DEBUG
    printf("\nMatrix[%d][%d]\n", n, m);
    printf("\n a string:\n");
    for(i=0; i<m-1; i++)
        printf("%d ",a[i]);
    printf("\n b string:\n");
    for(i=0; i<n-1; i++)
        printf("%d ",b[i]);
    printf("\n");
    #endif

    

    int Vsize = 8;

    double t;
    int it;
    for(it=0; it<NumOfTest; it++){
    long long int maxPos         = 0;
    long long int maxPos_max_len = 0;

    long long int ind   = 3;
    long long int indd  = 0;
    long long int indul = 1;
    long long int ind_u, ind_d, ind_l; 
   __m256i offset =_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
   __m256i reverseIndices = _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    for (i = 2; i < m+n-1; i++) { //Lines
        long long int max_len;
        long long int ii,jj;
        long long int j_start, j_end;
        if (i<n){
		    max_len = i+1;
            j_start = 1;
            j_end   = max_len-1;
            ind_u   = ind - max_len;
            ind_l   = ind - max_len + 1;
            ind_d   = ind - (max_len<<1) + 2;
            ii = i-j_start;
            jj = j_start;
        }
	    else if (i>=m){
		    max_len = m+n-1-i;
            j_start = 0;   
            j_end   = max_len;     
            ind_u   = ind - max_len - 1;
            ind_l   = ind - max_len; 
            ind_d   = ind - (max_len<<1) - 2;
            ii = m-1-j_start;
            jj = i-m+j_start+1;
        }
	    else{
		    max_len   = n;
            j_start = 1;
            j_end   = max_len;
            ind_u     = ind - max_len - 1;
            ind_l     = ind - max_len;
            ii = i-j_start;
            jj = j_start; 
            if(i>n)
                ind_d = ind - (max_len<<1) - 1;
            else
                ind_d = ind - (max_len<<1);
        }  
       
        __m256i* Hu = (__m256i*) (H+ind_u+j_start);
        __m256i* Hl = (__m256i*) (H+ind_l+j_start);
        __m256i* Hd = (__m256i*) (H+ind_d+j_start);
        __m256i* HH = (__m256i*) (H+ind+j_start);
        __m256i* PP = (__m256i*) (P+ind+j_start);
        
        for (j=j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns          
           similarityScoreIntrinsic(HH, Hu, Hd, Hl, PP, reverseIndices, ii, jj, H, ind+j, max_len, &maxPos, &maxPos_max_len, gapScore, a, b, m, n);
           ii -= Vsize;
           jj += Vsize;
           Hu++;
           Hl++;
           Hd++;
           HH++;
           PP++;
        }
        for(;j<j_end; j++){   
            similarityScore(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H, P, max_len, &maxPos, &maxPos_max_len, gapScore, a, b, m, n);
            ii --;
            jj ++;
        }
        ind += max_len;
    } 
       backtrack(P, maxPos, maxPos_max_len, m, n);
    }

    #ifdef DEBUG
    saveInFile(H, a, b, m, n);

    printf("\nSimilarity Matrix:\n");
    printMatrix(H, a, b, m, n);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P, a, b, m, n);
    #endif

}

void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, int* H, int* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, int gapScore, int8_t *a, int8_t *b, int m, int n) {

    int up, left, diag;

    //Get element above
    up = H[ind_u] + gapScore;

    //Get element on the left
    left = H[ind_l] + gapScore;

    //Get element on the diagonal
    diag = H[ind_d] + iBlosum62[a[ii-1]*32+b[jj-1]];

    //Calculates the maximum
    int max = NONE;
    int pred = NONE;
    /* === Matrix ===
     *      a[0] ... a[n] 
     * b[0]
     * ...
     * b[n]
     *
     * generate 'a' from 'b', if '←' insert e '↑' remove
     * a=GAATTCA
     * b=GACTT-A
     * 
     * generate 'b' from 'a', if '←' insert e '↑' remove
     * b=GACTT-A
     * a=GAATTCA
    */
    
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
    H[ind] = max;
    P[ind] = pred;

    //Updates maximum score to be used as seed on backtrack 
    if (max > H[*maxPos]) {
        *maxPos = ind;
        *maxPos_max_len = max_len;
    }

}  /* End of similarityScore */


void similarityScoreIntrinsic(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, long long int ii, long long int jj, int* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, int gapScore, int8_t *a, int8_t *b, int m, int n) {

   __m256i up, left, diag;

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);

    //Get element above
    up                    =_mm256_add_epi32(HHu,_mm256_set1_epi32(gapScore));

    //Get element on the left
    left                  =_mm256_add_epi32(HHl,_mm256_set1_epi32(gapScore));

    //Get element on the diagonal
    __m128i input = _mm_loadu_si128((__m128i*)(a+ii-8));
    __m256i A     = _mm256_cvtepu8_epi32(input);
            A     = _mm256_permutevar8x32_epi32(A, reverseIndices);
            input = _mm_loadu_si128((__m128i*)(b+jj-1));
    __m256i B     = _mm256_cvtepu8_epi32(input);

    __m256i addresses = _mm256_add_epi32(_mm256_mullo_epi32(A, _mm256_set1_epi32(32)), B);
    // Gather the values from the other matrix using the calculated addresses
    __m256i gatheredData = _mm256_i32gather_epi32((void*) iBlosum62, addresses, sizeof(int ));
     diag                  =_mm256_add_epi32(HHd, gatheredData);

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi32(NONE);
   __m256i pred =_mm256_set1_epi32(NONE);

    /* === Matrix ===
     *      a[0] ... a[n] 
     * b[0]
     * ...
     * b[n]
     *
     * generate 'a' from 'b', if '←' insert e '↑' remove
     * a=GAATTCA
     * b=GACTT-A
     * 
     * generate 'b' from 'a', if '←' insert e '↑' remove
     * b=GACTT-A
     * a=GAATTCA
    */
   //same letter ↖
    __m256i mask    = _mm256_cmpgt_epi32(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(DIAGONAL), mask);

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi32(up, max);
    max     = _mm256_blendv_epi8(max, up, mask);
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(UP), mask);

    //insert letter ←
    mask    = _mm256_cmpgt_epi32(left, max);
    max     = _mm256_blendv_epi8(max, left, mask);
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(LEFT), mask);

    //Inserts the value in the similarity and predecessor matrixes
    _mm256_storeu_si256(HH, max);
    _mm256_storeu_si256(PP, pred);
    
    //Updates maximum score to be used as seed on backtrack 
    __m256i vmax = max;
    vmax = _mm256_max_epu32(vmax, _mm256_alignr_epi8(vmax, vmax, 4));
    vmax = _mm256_max_epu32(vmax, _mm256_alignr_epi8(vmax, vmax, 8));
    vmax = _mm256_max_epu32(vmax, _mm256_permute2x128_si256(vmax, vmax, 0x01));

    __m256i vcmp = _mm256_cmpeq_epi32(max, vmax);

    int max_index = _mm256_movemask_epi8(vcmp);

    max_index = __builtin_ctz(max_index) >> 2;

    if (H[ind+max_index] > H[*maxPos]) {
        *maxPos         = ind+max_index;
        *maxPos_max_len = max_len;
    }

}  /* End of similarityScore */

/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(int* P, long long int maxPos, long long int maxPos_max_len, int m, int n) {
    //hold maxPos value
    long long int predPos;
    long long int predMaxLen;
    #ifdef pragmas
    #pragma GCC ivdep
    #endif
    //backtrack from maxPos to startPos = 0 
    long long int first_sec = (n*(n+1))/2;
    long long int last_sec  = n*m - (n*(n-1))/2;
    long long int ind_u, ind_d, ind_l; 
    bool diagCompensate     = 0;
    do {
        if (maxPos<first_sec){
            if(diagCompensate){
                if(maxPos<first_sec-n)
                    maxPos_max_len --;
                diagCompensate = 0;
            }
            ind_u      = maxPos - maxPos_max_len;
            ind_l      = maxPos - maxPos_max_len + 1;
            ind_d      = maxPos - (maxPos_max_len<<1) + 2;
            predMaxLen = maxPos_max_len-1;
        }
	    else if (maxPos>=last_sec){
            if(diagCompensate){
                maxPos_max_len ++;
                diagCompensate = 0;
            }
            ind_u      = maxPos - maxPos_max_len - 1;
            ind_l      = maxPos - maxPos_max_len; 
            ind_d      = maxPos - (maxPos_max_len<<1) - 2;
            predMaxLen = maxPos_max_len+1;
        }
	    else{
            if(diagCompensate){
                if(maxPos>=last_sec-n)
                    maxPos_max_len ++;
                diagCompensate = 0;
            }
            ind_u      = maxPos - n - 1;
            ind_l      = maxPos - n;
            predMaxLen = maxPos_max_len; 
            if(maxPos>=first_sec+n)
                ind_d  = maxPos - (n<<1) - 1;
            else
                ind_d  = maxPos - (n<<1);
        }

        if(P[maxPos] == DIAGONAL){
            predPos        = ind_d;
            diagCompensate = 1;
        }
        else if(P[maxPos] == UP)
            predPos        = ind_u;
        else if(P[maxPos] == LEFT)
            predPos        = ind_l;
        
        P[maxPos]*=PATH;
        maxPos         = predPos;
        maxPos_max_len = predMaxLen;
    } while(P[maxPos] != NONE);
}  /* End of backtrack */