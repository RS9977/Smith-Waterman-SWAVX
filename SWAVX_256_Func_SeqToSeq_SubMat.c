#include "SWAVX_256_Func_SeqToSeq_SubMat.h"

void SWAVX_256_SeqToSeq_SubMat(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal){

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

    #ifdef AFFINE
    INT* f = calloc(m, sizeof(INT));
    INT* e = calloc(n, sizeof(INT));
    #endif
    
    __m256i *local_max = malloc(sizeof(__m256i));
    
    #ifdef L8
    __m256i temp  = _mm256_set_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0);
    int Vsize = 32;
    __m256i reverseIndices2 = _mm256_set_epi8(0,   1,  2,  3,  4,  5,  6,  7,
                                             8,   9, 10, 11, 12, 13, 14, 15,
                                             16, 17, 18, 19, 20, 21, 22, 23,
                                             24, 25, 26, 27, 28, 29, 30, 31);
    #ifdef SUBMAT
    __m128i reverseIndices = _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    #else
    __m256i reverseIndices = _mm256_set_epi8(0,   1,  2,  3,  4,  5,  6,  7,
                                             8,   9, 10, 11, 12, 13, 14, 15,
                                             16, 17, 18, 19, 20, 21, 22, 23,
                                             24, 25, 26, 27, 28, 29, 30, 31);    
    #endif
    
    #elif L16
  
    __m256i temp  = _mm256_set_epi16(0, 0, 0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0, 0, 0);

    int Vsize = 16;
    __m128i reverseIndices = _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    __m256i reverseIndices2 = _mm256_set_epi8(1,   0,  3,  2,  5,  4,  7,  6,
                                              9,   8, 11, 10, 13, 12, 15, 14,
                                              17, 16, 19, 18, 21, 20, 23, 22,
                                              25, 24, 27, 26, 29, 28, 31, 30);
    #else
    //int* local_max_mem = malloc(8*sizeof(int));
    //__m256i *local_max = (__m256i*) local_max_mem;
    
    __m256i temp = _mm256_set_epi32(0, 0, 0, 0, 0, 0, 0, 0);
    
    int Vsize = 8;
    __m256i reverseIndices = _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    __m256i reverseIndices2 = _mm256_set_epi8(3, 2, 1, 0,  7, 6, 5, 4,
                                              11,   10, 9, 8, 15, 14, 13, 12,
                                              19, 18, 17, 16, 23, 22, 21, 20,
                                              27, 26, 25, 24, 31, 30, 29, 28);
    #endif
    _mm256_storeu_si256(local_max, temp);

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
            #ifndef SAVEHP
            *(H+ind) = 0;
            *(H+ind+max_len-1) = 0;
            #endif
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
            #ifndef SAVEHP
            *(H+ind) = 0;
            #endif
            if(i>n)
                ind_d = ind - (max_len<<1) - 1;
            else
                ind_d = ind - (max_len<<1);
        }  
      
        __m256i* Hu = (__m256i*) (H+ind_u+j_start);
        __m256i* Hl = (__m256i*) (H+ind_l+j_start);
        __m256i* Hd = (__m256i*) (H+ind_d+j_start);
        __m256i* HH = (__m256i*) (H+ind+j_start);
        __m256i* PP;
        #ifdef BT
        PP = (__m256i*) (P+ind+j_start);
        #endif

        #ifndef AFFINE
        for (j=j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns          
           
           #ifdef L8
           similarityScoreIntrinsic8(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, m, n, 32);
           #elif L16
           similarityScoreIntrinsic16(HH, Hu, Hd, Hl, PP, reverseIndices, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, m, n);
           #else
           similarityScoreIntrinsic32(HH, Hu, Hd, Hl, PP, reverseIndices, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, m, n);
           #endif

           ii -= Vsize;
           jj += Vsize;
           Hu++;
           Hl++;
           Hd++;
           HH++;
           PP++;
        }
        #ifdef L8
            
            if(j_end-j > 4){
                similarityScoreIntrinsic8(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, m, n, j_end-j-1);
                j = j_end;
            }
            
        #endif
        for(;j<j_end; j++){   
            similarityScore(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H, P, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, m, n);
            ii --;
            jj ++;
        }
        #else
        for (j=j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns          
           
           #ifdef L8
           similarityScoreIntrinsic8_affine(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, e, f, m, n, 32);
           #elif L16
           similarityScoreIntrinsic16_affine(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, e, f, m, n, 16);
           #else
           similarityScoreIntrinsic32_affine(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, e, f, m, n);
           #endif

           ii -= Vsize;
           jj += Vsize;
           Hu++;
           Hl++;
           Hd++;
           HH++;
           PP++;
        }
        #ifdef L8
        if(j_end-j > 4){
            similarityScoreIntrinsic8_affine(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, e, f, m, n, j_end-j-1);
            j = j_end;
        }
        #elif L16
        if(j_end-j > 4){
            similarityScoreIntrinsic16_affine(HH, Hu, Hd, Hl, PP, reverseIndices, reverseIndices2, ii, jj, H, ind+j, local_max, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, e, f, m, n, j_end-j-1);
            j = j_end;
        }
        #endif
        for(;j<j_end; j++){   
            similarityScore_affine(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H, P, max_len, &maxPos, &maxPos_max_len, maxVal, a, b, e, f, m, n);
            ii --;
            jj ++;
        }
        #endif
        ind += max_len;
    } 
    
       #ifdef BT
       backtrack(P, maxPos, maxPos_max_len, m, n);
       #endif
    }

    #ifndef BT
    __m256i max = _mm256_loadu_si256(local_max);
    #ifdef L8
    __m256i vtmp1;
    __m256i vtmp2;
    vtmp1 = _mm256_permute2x128_si256(max, max, 1);
    vtmp1 = _mm256_max_epi8(vtmp1, max);
    vtmp2 = _mm256_srli_si256(vtmp1, 8);
    vtmp1 = _mm256_max_epi8(vtmp1, vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,4);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,2);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,1);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    short int v = _mm256_extract_epi8(vtmp1,0);
    *maxVal = (v>*maxVal)? v: *maxVal;
    #elif L16
    __m256i vtmp1;
    __m256i vtmp2;
    
    vtmp1 = _mm256_permute2x128_si256(max, max, 1);
    vtmp1 = _mm256_max_epi16(vtmp1, max);
    vtmp2 =  _mm256_permute4x64_epi64(vtmp1, 0x01);
    vtmp1 = _mm256_max_epi16(vtmp1, vtmp2);
    vtmp2 = _mm256_permutevar8x32_epi32(vtmp1, _mm256_setr_epi32(1, 0, 0, 0, 0, 0, 0, 0));
    vtmp1 = _mm256_max_epi16(vtmp1, vtmp2);
    short int v1 = _mm256_extract_epi16(vtmp1,0);
    short int v2 = _mm256_extract_epi16(vtmp1,1);
    v1 = (v1>v2)?v1:v2;
    *maxVal = (v1>*maxVal)? v1: *maxVal;
    
    #else 
    __m256i vmax = max;
    vmax = _mm256_max_epu32(vmax, _mm256_alignr_epi8(vmax, vmax, 4));
    vmax = _mm256_max_epu32(vmax, _mm256_alignr_epi8(vmax, vmax, 8));
    vmax = _mm256_max_epu32(vmax, _mm256_permute2x128_si256(vmax, vmax, 0x01));

    __m256i vcmp = _mm256_cmpeq_epi32(max, vmax);

    int max_index = _mm256_movemask_epi8(vcmp);

    max_index = __builtin_ctz(max_index) >> 2;
    int v = _mm256_extract_epi32(vmax, 0);
    *maxVal = (v>*maxVal)? v: *maxVal;

    #endif
    #endif

    #ifdef DEBUG
    saveInFile(H, a, b, m, n);

    printf("\nSimilarity Matrix:\n");
    printMatrix(H, a, b, m, n);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P, a, b, m, n);
    #endif

    #ifdef AFFINE
    free(f);
    free(e);
    #endif

    free(local_max);
}

void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, INT* H, INT* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, int m, int n) {

    INT up, left, diag;

    //Get element above
    up   = H[ind_u] + gapScore;

    //Get element on the left
    left = H[ind_l] + gapScore;

    #ifdef SUBMAT
    //Get element on the diagonal
    diag = H[ind_d] + iBlosum62[a[ii-1]*32+b[jj-1]];
    #else
    diag = H[ind_d] + matchMissmatchScore(ii, jj, a, b);
    #endif

    //Calculates the maximum
    INT max = NONE;
    #ifdef BT
    INT pred = NONE;
    #endif
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
        #ifdef BT
        pred = DIAGONAL;
        #endif
    }

    if (up > max) { //remove letter ↑ 
        max = up;
        #ifdef BT
        pred = UP;
        #endif
    }
    
    if (left > max) { //insert letter ←
        max = left;
        #ifdef BT
        pred = LEFT;
        #endif
    }
    //Inserts the value in the similarity and predecessor matrixes
    H[ind] = max;
    #ifdef BT
    P[ind] = pred;
    #endif

    //Updates maximum score to be used as seed on backtrack 
    #ifdef BT
    if (max > H[*maxPos]) {
        *maxPos = ind;
        *maxPos_max_len = max_len;
    }
    #else
    *maxVal = (max>*maxVal)? max: *maxVal;
    #endif

}  /* End of similarityScore */


void similarityScoreIntrinsic32(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, long long int ii, long long int jj, INT* H, long long int ind, __m256i *local_max  , long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, int m, int n) {

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

    __m256i mask; 
    #ifdef SUBMAT
    __m256i addresses = _mm256_add_epi32(_mm256_mullo_epi32(A, _mm256_set1_epi32(32)), B);
    // Gather the values from the other matrix using the calculated addresses
    __m256i gatheredData   = _mm256_i32gather_epi32((void*) iBlosum62, addresses, sizeof(int ));
     diag                  = _mm256_add_epi32(HHd, gatheredData);

    #else
            mask           = _mm256_cmpeq_epi32(A, B);

    __m256i MATCHSCORE     = _mm256_set1_epi32(matchScore);
    __m256i MISSMATCHSCORE = _mm256_set1_epi32(missmatchScore);
    __m256i MATCHMISS      = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                   = _mm256_add_epi32(HHd, MATCHMISS);
    #endif

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi32(NONE);
   #ifdef BT
   __m256i pred =_mm256_set1_epi32(NONE);
   #endif


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
    mask    = _mm256_cmpgt_epi32(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(DIAGONAL), mask);
    #endif

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi32(up, max);
    max     = _mm256_blendv_epi8(max, up, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(UP), mask);
    #endif

    //insert letter ←
    mask    = _mm256_cmpgt_epi32(left, max);
    max     = _mm256_blendv_epi8(max, left, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(LEFT), mask);
    #endif

    //Inserts the value in the similarity and predecessor matrixes
    _mm256_storeu_si256(HH, max);
    #ifdef BT
    _mm256_storeu_si256(PP, pred);
    #endif
    
    //Updates maximum score to be used as seed on backtrack 
    #ifdef BT
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
    #else
    //int v = _mm256_extract_epi32(vmax, 0);
    //*maxVal = (v>*maxVal)? v: *maxVal;
    __m256i tempmax = _mm256_loadu_si256(local_max);
    tempmax = _mm256_max_epi32(max, tempmax);
    _mm256_storeu_si256(local_max,tempmax);
    #endif

}  /* End of similarityScore */


void similarityScoreIntrinsic16(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m128i reverseIndices, long long int ii, long long int jj, INT* H, long long int ind, __m256i *local_max   , long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, int m, int n) {

   __m256i up, left, diag;

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);

    //Get element above
    up                    =_mm256_add_epi16(HHu,_mm256_set1_epi16(gapScore));

    //Get element on the left
    left                  =_mm256_add_epi16(HHl,_mm256_set1_epi16(gapScore));

    //Get element on the diagonal

    __m128i input = _mm_loadu_si128((__m128i*)(a+ii-16));
            input = _mm_shuffle_epi8(input, reverseIndices);
    __m256i A     = _mm256_cvtepu8_epi16(input);
            input = _mm_loadu_si128((__m128i*)(b+jj-1));
    __m256i B     = _mm256_cvtepu8_epi16(input);    

    __m256i mask; 
    #ifdef SUBMAT
    __m256i addresses        = _mm256_add_epi16(_mm256_mullo_epi16(A, _mm256_set1_epi16(32)), B);
    __m256i lowBits          = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 0));
    __m256i highBits         = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 1)); 
    __m256i gatheredDataLow  = _mm256_i32gather_epi32((void*) iBlosum62, lowBits,  sizeof(int));
    __m256i gatheredDataHigh = _mm256_i32gather_epi32((void*) iBlosum62, highBits, sizeof(int));
    __m256i gatheredData     = _mm256_packs_epi32(gatheredDataLow, gatheredDataHigh);
    gatheredData             = _mm256_permute4x64_epi64(gatheredData, 0xd8);
     diag                    = _mm256_add_epi16(HHd, gatheredData);
    #else
            mask           = _mm256_cmpeq_epi16(A, B);    
    __m256i MATCHSCORE     = _mm256_set1_epi16(matchScore);
    __m256i MISSMATCHSCORE = _mm256_set1_epi16(missmatchScore);
    __m256i MATCHMISS      = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                   = _mm256_add_epi16(HHd, MATCHMISS);
    #endif

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi16(NONE);
   #ifdef BT
   __m256i pred =_mm256_set1_epi16(NONE);
   #endif


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
    mask    = _mm256_cmpgt_epi16(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi16(DIAGONAL), mask);
    #endif

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi16(up, max);
    max     = _mm256_blendv_epi8(max, up, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi16(UP), mask);
    #endif

    //insert letter ←
    mask    = _mm256_cmpgt_epi16(left, max);
    max     = _mm256_blendv_epi8(max, left, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi16(LEFT), mask);
    #endif

    //Inserts the value in the similarity and predecessor matrixes
    _mm256_storeu_si256(HH, max);
    #ifdef BT
    _mm256_storeu_si256(PP, pred);
    #endif
    

    #ifdef BT
    __m256i vtmp1;
    __m256i vtmp2;
    
    vtmp1 = _mm256_permute2x128_si256(max, max, 1);
    vtmp1 = _mm256_max_epi16(vtmp1, max);
    vtmp2 =  _mm256_permute4x64_epi64(vtmp1, 0x01);
    vtmp1 = _mm256_max_epi16(vtmp1, vtmp2);
    vtmp2 = _mm256_permutevar8x32_epi32(vtmp1, _mm256_setr_epi32(1, 0, 0, 0, 0, 0, 0, 0));
    vtmp1 = _mm256_max_epi16(vtmp1, vtmp2);
    short int v1 = _mm256_extract_epi16(vtmp1,0);
    short int v2 = _mm256_extract_epi16(vtmp1,1);
    v1 = (v1>v2)?v1:v2;

    
    __m256i vcmp = _mm256_cmpeq_epi16(_mm256_set1_epi16(v1), max);
    int max_index = _mm256_movemask_epi8(vcmp);
    max_index = __builtin_ctz(max_index) >> 1;



    if (H[ind+max_index] > H[*maxPos]) {
        *maxPos         = ind+max_index;
        *maxPos_max_len = max_len;
    }
    #else
    //*maxVal = (v1>*maxVal)? v1: *maxVal;
    __m256i tempmax = _mm256_loadu_si256(local_max);
    tempmax = _mm256_max_epi16(max, tempmax);
    _mm256_storeu_si256(local_max,tempmax);
    #endif
}  /* End of similarityScore */



void similarityScoreIntrinsic8(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, RevType reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, __m256i *local_max   , long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, int m, int n, int k) {

   __m256i up, left, diag;
    

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);

    //Get element above
    up                    =_mm256_add_epi8(HHu,_mm256_set1_epi8(gapScore));

    //Get element on the left
    left                  =_mm256_add_epi8(HHl,_mm256_set1_epi8(gapScore));

    //Get element on the diagonal

    

    __m256i mask; 
    #ifdef SUBMAT
    __m128i input = _mm_loadu_si128((__m128i*)(a+ii-16));
            input = _mm_shuffle_epi8(input, reverseIndices);
    __m256i A     = _mm256_cvtepu8_epi16(input);
            input = _mm_loadu_si128((__m128i*)(b+jj-1));
    __m256i B     = _mm256_cvtepu8_epi16(input);
        
    __m256i addresses        = _mm256_add_epi16(_mm256_mullo_epi16(A, _mm256_set1_epi16(32)), B);
    if(k<16){
        __m256i masktemp = _mm256_cmpgt_epi16(_mm256_set_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), _mm256_set1_epi16(16-k-1));
        addresses = _mm256_and_si256(addresses, masktemp);
    }
    __m256i lowBits          = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 0));
    __m256i highBits         = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 1)); 
    __m256i gatheredDataLow  = _mm256_i32gather_epi32((void*) iBlosum62, lowBits,  sizeof(int));
    __m256i gatheredDataHigh = _mm256_i32gather_epi32((void*) iBlosum62, highBits, sizeof(int));
    __m256i gatheredData1     = _mm256_packs_epi32(gatheredDataLow, gatheredDataHigh);
    gatheredData1             = _mm256_permute4x64_epi64(gatheredData1, 0xd8);
    input                     = _mm_loadu_si128((__m128i*)(a+ii-32));
            input = _mm_shuffle_epi8(input, reverseIndices);
    __m256i gatheredData2;
    if(k>=16){
     A     = _mm256_cvtepu8_epi16(input);
            input = _mm_loadu_si128((__m128i*)(b+jj+15));
     B     = _mm256_cvtepu8_epi16(input);    
     addresses        = _mm256_add_epi16(_mm256_mullo_epi16(A, _mm256_set1_epi16(32)), B);
     if(k<32){
        __m256i masktemp = _mm256_cmpgt_epi16(_mm256_set_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), _mm256_set1_epi16(32-k-1));
        addresses = _mm256_and_si256(addresses, masktemp);
    }
     lowBits          = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 0));
     highBits         = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 1)); 
     gatheredDataLow  = _mm256_i32gather_epi32((void*) iBlosum62, lowBits,  sizeof(int));
     gatheredDataHigh = _mm256_i32gather_epi32((void*) iBlosum62, highBits, sizeof(int));
    gatheredData2     = _mm256_packs_epi32(gatheredDataLow, gatheredDataHigh);
    gatheredData2             = _mm256_permute4x64_epi64(gatheredData2, 0xd8); 
    }
    else{
    gatheredData2 = _mm256_setzero_si256();
    }
    __m256i gatheredData     = _mm256_packs_epi16(gatheredData1,gatheredData2);
     gatheredData             = _mm256_permute4x64_epi64(gatheredData, 0xd8);
    diag                    = _mm256_add_epi8(HHd, gatheredData);
    #else
    __m256i input = _mm256_loadu_si256((__m256i*)(a+ii-32));
    __m256i A     = _mm256_shuffle_epi8(input, reverseIndices);
            A     = _mm256_permute2x128_si256(A,A,3);
    __m256i B     = _mm256_loadu_si256((__m256i*)(b+jj-1));   
            mask           = _mm256_cmpeq_epi8(A, B);    
    __m256i MATCHSCORE     = _mm256_set1_epi8(matchScore);
    __m256i MISSMATCHSCORE = _mm256_set1_epi8(missmatchScore);
    __m256i MATCHMISS      = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                   = _mm256_add_epi8(HHd, MATCHMISS);
    #endif

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi8(NONE);
   #ifdef BT
   __m256i pred =_mm256_set1_epi8(NONE);
   #endif


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
    mask    = _mm256_cmpgt_epi8(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi8(DIAGONAL), mask);
    #endif

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi8(up, max);
    max     = _mm256_blendv_epi8(max, up, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi8(UP), mask);
    #endif

    //insert letter ←
    mask    = _mm256_cmpgt_epi8(left, max);
    max     = _mm256_blendv_epi8(max, left, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi8(LEFT), mask);
    #endif

    //Inserts the value in the similarity and predecessor matrixes
    if(k<32){
        __m256i masktemp = _mm256_cmpgt_epi8(reverseIndices2, _mm256_set1_epi8(32-k-1));
        max = _mm256_and_si256(max, masktemp);
    }
    _mm256_storeu_si256(HH, max);
    #ifdef BT
    _mm256_storeu_si256(PP, pred);
    #endif
    
    #ifdef BT
    __m256i vtmp1;
    __m256i vtmp2;
    vtmp1 = _mm256_permute2x128_si256(max, max, 1);
    vtmp1 = _mm256_max_epi8(vtmp1, max);
    vtmp2 = _mm256_srli_si256(vtmp1, 8);
    vtmp1 = _mm256_max_epi8(vtmp1, vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,4);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,2);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,1);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    short int v = _mm256_extract_epi8(vtmp1,0);
    
    __m256i vcmp = _mm256_cmpeq_epi8(_mm256_set1_epi8(v), max);
    int max_index = _mm256_movemask_epi8(vcmp);
    max_index = __builtin_ctz(max_index);

    
    if (H[ind+max_index] > H[*maxPos]) {
        *maxPos         = ind+max_index;
        *maxPos_max_len = max_len;
    }
    #else
    //*maxVal = (v>*maxVal)? v: *maxVal;
    __m256i tempmax = _mm256_loadu_si256(local_max);
    tempmax = _mm256_max_epi8(max, tempmax);
    _mm256_storeu_si256(local_max,tempmax);
    #endif
}  /* End of similarityScore */


void similarityScore_affine(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, INT* H, INT* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n) {

    INT uph, upf, lefth, lefte, diag;

    //Get element above


    #ifdef SUBMAT
    //Get element on the diagonal
    diag = H[ind_d] + iBlosum62[a[ii-1]*32+b[jj-1]];
    #else
    diag = H[ind_d] + matchMissmatchScore(ii, jj, a, b);
    #endif

    //Calculates the maximum
    INT max = NONE;
    #ifdef BT
    INT pred = NONE;
    #endif
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
        #ifdef BT
        pred = DIAGONAL;
        #endif
    }

    if (f[ii-1] > max) { //remove letter ↑ 
        max = f[ii-1];
        #ifdef BT
        pred = UP;
        #endif
    }
    
    if (e[jj-1] > max) { //insert letter ←
        max = e[jj-1];
        #ifdef BT
        pred = LEFT;
        #endif
    }
    //Inserts the value in the similarity and predecessor matrixes
    H[ind] = max;
    f[ii-1] = (f[ii-1]+gapExt>max+gapOpen)?f[ii-1]+gapExt:max+gapOpen;
    e[jj-1] = (e[jj-1]+gapExt>max+gapOpen)?e[jj-1]+gapExt:max+gapOpen;
    #ifdef BT
    P[ind] = pred;
    #endif

    //Updates maximum score to be used as seed on backtrack 
    #ifdef BT
    if (max > H[*maxPos]) {
        *maxPos = ind;
        *maxPos_max_len = max_len;
    }
    #else
    *maxVal = (max>*maxVal)? max: *maxVal;
    #endif

}  /* End of similarityScore */

void similarityScoreIntrinsic32_affine(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, __m256i *local_max  , long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n) {

    __m256i upf, lefte, uph, lefth, diag;
    
    __m256i F = _mm256_loadu_si256((__m256i*)(f+ii-8));
            F = _mm256_shuffle_epi8(F, reverseIndices2);
            F = _mm256_permute2x128_si256(F,F,3);
    __m256i E = _mm256_loadu_si256((__m256i*)(e+jj-1));

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);
    

    //Get element on the diagonal
    __m128i input = _mm_loadu_si128((__m128i*)(a+ii-8));
    __m256i A     = _mm256_cvtepu8_epi32(input);
            A     = _mm256_permutevar8x32_epi32(A, reverseIndices);
            input = _mm_loadu_si128((__m128i*)(b+jj-1));
    __m256i B     = _mm256_cvtepu8_epi32(input);

    __m256i mask; 
    #ifdef SUBMAT
    __m256i addresses = _mm256_add_epi32(_mm256_mullo_epi32(A, _mm256_set1_epi32(32)), B);
    // Gather the values from the other matrix using the calculated addresses
    __m256i gatheredData   = _mm256_i32gather_epi32((void*) iBlosum62, addresses, sizeof(int ));
     diag                  = _mm256_add_epi32(HHd, gatheredData);

    #else
            mask           = _mm256_cmpeq_epi32(A, B);

    __m256i MATCHSCORE     = _mm256_set1_epi32(matchScore);
    __m256i MISSMATCHSCORE = _mm256_set1_epi32(missmatchScore);
    __m256i MATCHMISS      = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                   = _mm256_add_epi32(HHd, MATCHMISS);
    #endif

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi32(NONE);
   #ifdef BT
   __m256i pred =_mm256_set1_epi32(NONE);
   #endif


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
    mask    = _mm256_cmpgt_epi32(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(DIAGONAL), mask);
    #endif

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi32(F, max);
    max     = _mm256_blendv_epi8(max, F, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(UP), mask);
    #endif

    //insert letter ←
    mask    = _mm256_cmpgt_epi32(E, max);
    max     = _mm256_blendv_epi8(max, E, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi32(LEFT), mask);
    #endif

    uph                    =_mm256_add_epi32(max,_mm256_set1_epi32(gapOpen));
    lefth                  =_mm256_add_epi32(max,_mm256_set1_epi32(gapOpen));
    upf                    =_mm256_add_epi32(F,_mm256_set1_epi32(gapExt));
    lefte                  =_mm256_add_epi32(E,_mm256_set1_epi32(gapExt));
    F     = _mm256_shuffle_epi8(F, reverseIndices2);
    F     = _mm256_permute2x128_si256(F,F,3);
    _mm256_storeu_si256((__m256i*)(f+ii-8), F);
    _mm256_storeu_si256((__m256i*)(e+jj-1), E);

    //Inserts the value in the similarity and predecessor matrixes
    _mm256_storeu_si256(HH, max);
    #ifdef BT
    _mm256_storeu_si256(PP, pred);
    #endif
    
    //Updates maximum score to be used as seed on backtrack 
    #ifdef BT
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
    #else
    //int v = _mm256_extract_epi32(vmax, 0);
    //*maxVal = (v>*maxVal)? v: *maxVal;
    __m256i tempmax = _mm256_loadu_si256(local_max);
    tempmax = _mm256_max_epi32(max, tempmax);
    _mm256_storeu_si256(local_max,tempmax);
    #endif

}  /* End of similarityScore */


void similarityScoreIntrinsic16_affine(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m128i reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, __m256i *local_max   , long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n, int k) {

    __m256i upf, lefte, uph, lefth, diag;
    

    __m256i F = _mm256_loadu_si256((__m256i*)(f+ii-16));
            F = _mm256_shuffle_epi8(F, reverseIndices2);
            F = _mm256_permute2x128_si256(F,F,3);
    __m256i E = _mm256_loadu_si256((__m256i*)(e+jj-1));

    

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);

   
    __m128i input = _mm_loadu_si128((__m128i*)(a+ii-16));
            input = _mm_shuffle_epi8(input, reverseIndices);
    __m256i A     = _mm256_cvtepu8_epi16(input);
            input = _mm_loadu_si128((__m128i*)(b+jj-1));
    __m256i B     = _mm256_cvtepu8_epi16(input);    

    __m256i mask; 
    #ifdef SUBMAT
    __m256i addresses        = _mm256_add_epi16(_mm256_mullo_epi16(A, _mm256_set1_epi16(32)), B);
    if(k<16){
    __m256i masktemp         = _mm256_cmpgt_epi16(_mm256_set_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), _mm256_set1_epi16(16-k-1));
            addresses        = _mm256_and_si256(addresses, masktemp);
    }
    __m256i lowBits          = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 0));
    __m256i highBits         = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 1)); 
    __m256i gatheredDataLow  = _mm256_i32gather_epi32((void*) iBlosum62, lowBits,  sizeof(int));
    __m256i gatheredDataHigh = _mm256_i32gather_epi32((void*) iBlosum62, highBits, sizeof(int));
    __m256i gatheredData     = _mm256_packs_epi32(gatheredDataLow, gatheredDataHigh);
    gatheredData             = _mm256_permute4x64_epi64(gatheredData, 0xd8);
     diag                    = _mm256_add_epi16(HHd, gatheredData);
    #else
            mask           = _mm256_cmpeq_epi16(A, B);    
    __m256i MATCHSCORE     = _mm256_set1_epi16(matchScore);
    __m256i MISSMATCHSCORE = _mm256_set1_epi16(missmatchScore);
    __m256i MATCHMISS      = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                   = _mm256_add_epi16(HHd, MATCHMISS);
    #endif

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi16(NONE);
   #ifdef BT
   __m256i pred =_mm256_set1_epi16(NONE);
   #endif


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
    mask    = _mm256_cmpgt_epi16(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi16(DIAGONAL), mask);
    #endif

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi16(F, max);
    max     = _mm256_blendv_epi8(max, F, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi16(UP), mask);
    #endif

    //insert letter ←
    mask    = _mm256_cmpgt_epi16(E, max);
    max     = _mm256_blendv_epi8(max, E, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi16(LEFT), mask);
    #endif

    uph                    =_mm256_add_epi16(max,_mm256_set1_epi16(gapOpen));
    lefth                  =_mm256_add_epi16(max,_mm256_set1_epi16(gapOpen));
    upf                    =_mm256_add_epi16(F,_mm256_set1_epi16(gapExt));
    lefte                  =_mm256_add_epi16(E,_mm256_set1_epi16(gapExt));
    F     = _mm256_shuffle_epi8(F, reverseIndices2);
    F     = _mm256_permute2x128_si256(F,F,3);
    _mm256_storeu_si256((__m256i*)(f+ii-16), F);
    _mm256_storeu_si256((__m256i*)(e+jj-1), E);
    //Inserts the value in the similarity and predecessor matrixes
    if(k<16){
        __m256i masktemp = _mm256_cmpgt_epi16(reverseIndices2, _mm256_set1_epi16(16-k-1));
        max = _mm256_and_si256(max, masktemp);
    }

    //Inserts the value in the similarity and predecessor matrixes
    _mm256_storeu_si256(HH, max);
    #ifdef BT
    _mm256_storeu_si256(PP, pred);
    #endif
    

    #ifdef BT
    __m256i vtmp1;
    __m256i vtmp2;
    
    vtmp1 = _mm256_permute2x128_si256(max, max, 1);
    vtmp1 = _mm256_max_epi16(vtmp1, max);
    vtmp2 =  _mm256_permute4x64_epi64(vtmp1, 0x01);
    vtmp1 = _mm256_max_epi16(vtmp1, vtmp2);
    vtmp2 = _mm256_permutevar8x32_epi32(vtmp1, _mm256_setr_epi32(1, 0, 0, 0, 0, 0, 0, 0));
    vtmp1 = _mm256_max_epi16(vtmp1, vtmp2);
    short int v1 = _mm256_extract_epi16(vtmp1,0);
    short int v2 = _mm256_extract_epi16(vtmp1,1);
    v1 = (v1>v2)?v1:v2;

    
    __m256i vcmp = _mm256_cmpeq_epi16(_mm256_set1_epi16(v1), max);
    int max_index = _mm256_movemask_epi8(vcmp);
    max_index = __builtin_ctz(max_index) >> 1;



    if (H[ind+max_index] > H[*maxPos]) {
        *maxPos         = ind+max_index;
        *maxPos_max_len = max_len;
    }
    #else
    //*maxVal = (v1>*maxVal)? v1: *maxVal;
    __m256i tempmax = _mm256_loadu_si256(local_max);
    tempmax = _mm256_max_epi16(max, tempmax);
    _mm256_storeu_si256(local_max,tempmax);
    #endif
}  /* End of similarityScore */


void similarityScoreIntrinsic8_affine(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, RevType reverseIndices, __m256i reverseIndices2, long long int ii, long long int jj, INT* H, long long int ind, __m256i *local_max   , long long int max_len, long long int* maxPos, long long int *maxPos_max_len, INT* maxVal, int8_t *a, int8_t *b, INT* e, INT* f, int m, int n, int k) {

   __m256i upf, lefte, uph, lefth, diag;
    
    __m256i F  = _mm256_loadu_si256((__m256i*)(f+ii-32));
            F     = _mm256_shuffle_epi8(F, reverseIndices2);
            F     = _mm256_permute2x128_si256(F,F,3);
    __m256i E  = _mm256_loadu_si256((__m256i*)(e+jj-1));

    __m256i HHu = _mm256_loadu_si256(Hu);
    __m256i HHd = _mm256_loadu_si256(Hd);
    __m256i HHl = _mm256_loadu_si256(Hl);

    //Get element on the diagonal

    

    __m256i mask; 
    #ifdef SUBMAT
    __m128i input = _mm_loadu_si128((__m128i*)(a+ii-16));
            input = _mm_shuffle_epi8(input, reverseIndices);
    __m256i A     = _mm256_cvtepu8_epi16(input);
            input = _mm_loadu_si128((__m128i*)(b+jj-1));
    __m256i B     = _mm256_cvtepu8_epi16(input);
        
    __m256i addresses        = _mm256_add_epi16(_mm256_mullo_epi16(A, _mm256_set1_epi16(32)), B);
    if(k<16){
        __m256i masktemp = _mm256_cmpgt_epi16(_mm256_set_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), _mm256_set1_epi16(16-k-1));
        addresses = _mm256_and_si256(addresses, masktemp);
    }
    __m256i lowBits          = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 0));
    __m256i highBits         = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 1)); 
    __m256i gatheredDataLow  = _mm256_i32gather_epi32((void*) iBlosum62, lowBits,  sizeof(int));
    __m256i gatheredDataHigh = _mm256_i32gather_epi32((void*) iBlosum62, highBits, sizeof(int));
    __m256i gatheredData1     = _mm256_packs_epi32(gatheredDataLow, gatheredDataHigh);
    gatheredData1             = _mm256_permute4x64_epi64(gatheredData1, 0xd8);
    input                     = _mm_loadu_si128((__m128i*)(a+ii-32));
            input = _mm_shuffle_epi8(input, reverseIndices);
    __m256i gatheredData2;
    if(k>=16){
     A     = _mm256_cvtepu8_epi16(input);
            input = _mm_loadu_si128((__m128i*)(b+jj+15));
     B     = _mm256_cvtepu8_epi16(input);    
     addresses        = _mm256_add_epi16(_mm256_mullo_epi16(A, _mm256_set1_epi16(32)), B);
     if(k<32){
        __m256i masktemp = _mm256_cmpgt_epi16(_mm256_set_epi16(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15), _mm256_set1_epi16(32-k-1));
        addresses = _mm256_and_si256(addresses, masktemp);
    }
     lowBits          = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 0));
     highBits         = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(addresses, 1)); 
     gatheredDataLow  = _mm256_i32gather_epi32((void*) iBlosum62, lowBits,  sizeof(int));
     gatheredDataHigh = _mm256_i32gather_epi32((void*) iBlosum62, highBits, sizeof(int));
    gatheredData2     = _mm256_packs_epi32(gatheredDataLow, gatheredDataHigh);
    gatheredData2             = _mm256_permute4x64_epi64(gatheredData2, 0xd8); 
    }
    else{
    gatheredData2 = _mm256_setzero_si256();
    }
    __m256i gatheredData     = _mm256_packs_epi16(gatheredData1,gatheredData2);
     gatheredData             = _mm256_permute4x64_epi64(gatheredData, 0xd8);
    diag                    = _mm256_add_epi8(HHd, gatheredData);
    #else
    __m256i input = _mm256_loadu_si256((__m256i*)(a+ii-32));
    __m256i A     = _mm256_shuffle_epi8(input, reverseIndices);
            A     = _mm256_permute2x128_si256(A,A,3);
    __m256i B     = _mm256_loadu_si256((__m256i*)(b+jj-1));   
            mask           = _mm256_cmpeq_epi8(A, B);    
    __m256i MATCHSCORE     = _mm256_set1_epi8(matchScore);
    __m256i MISSMATCHSCORE = _mm256_set1_epi8(missmatchScore);
    __m256i MATCHMISS      = _mm256_blendv_epi8(MISSMATCHSCORE, MATCHSCORE, mask);
    diag                   = _mm256_add_epi8(HHd, MATCHMISS);
    #endif

    //Calculates the maximum
   __m256i max  =_mm256_set1_epi8(NONE);
   #ifdef BT
   __m256i pred =_mm256_set1_epi8(NONE);
   #endif


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
    mask    = _mm256_cmpgt_epi8(diag, max);
    max     = _mm256_blendv_epi8(max, diag, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi8(DIAGONAL), mask);
    #endif

    //remove letter ↑ 
    mask    = _mm256_cmpgt_epi8(F, max);
    max     = _mm256_blendv_epi8(max, F, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi8(UP), mask);
    #endif

    //insert letter ←
    mask    = _mm256_cmpgt_epi8(E, max);
    max     = _mm256_blendv_epi8(max, E, mask);
    #ifdef BT
    pred    = _mm256_blendv_epi8(pred, _mm256_set1_epi8(LEFT), mask);
    #endif

    uph                    =_mm256_add_epi8(max,_mm256_set1_epi8(gapOpen));
    lefth                  =_mm256_add_epi8(max,_mm256_set1_epi8(gapOpen));
    upf                    =_mm256_add_epi8(F,_mm256_set1_epi8(gapExt));
    lefte                  =_mm256_add_epi8(E,_mm256_set1_epi8(gapExt));
    F     = _mm256_shuffle_epi8(F, reverseIndices2);
    F     = _mm256_permute2x128_si256(F,F,3);
    _mm256_storeu_si256((__m256i*)(f+ii-32), F);
    _mm256_storeu_si256((__m256i*)(e+jj-1), E);
    //Inserts the value in the similarity and predecessor matrixes
    if(k<32){
        __m256i masktemp = _mm256_cmpgt_epi8(reverseIndices2, _mm256_set1_epi8(32-k-1));
        max = _mm256_and_si256(max, masktemp);
    }

    _mm256_storeu_si256(HH, max);
    #ifdef BT
    _mm256_storeu_si256(PP, pred);
    #endif
    
    #ifdef BT
    __m256i vtmp1;
    __m256i vtmp2;
    vtmp1 = _mm256_permute2x128_si256(max, max, 1);
    vtmp1 = _mm256_max_epi8(vtmp1, max);
    vtmp2 = _mm256_srli_si256(vtmp1, 8);
    vtmp1 = _mm256_max_epi8(vtmp1, vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,4);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,2);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    vtmp2 = _mm256_srli_si256(vtmp1,1);
    vtmp1 = _mm256_max_epi8(vtmp1,vtmp2);
    short int v = _mm256_extract_epi8(vtmp1,0);
    
    __m256i vcmp = _mm256_cmpeq_epi8(_mm256_set1_epi8(v), max);
    int max_index = _mm256_movemask_epi8(vcmp);
    max_index = __builtin_ctz(max_index);

    
    if (H[ind+max_index] > H[*maxPos]) {
        *maxPos         = ind+max_index;
        *maxPos_max_len = max_len;
    }
    #else
    //*maxVal = (v>*maxVal)? v: *maxVal;
    __m256i tempmax = _mm256_loadu_si256(local_max);
    tempmax = _mm256_max_epi8(max, tempmax);
    _mm256_storeu_si256(local_max,tempmax);
    #endif
}  /* End of similarityScore */

/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(INT* P, long long int maxPos, long long int maxPos_max_len, int m, int n) {
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

/*--------------------------------------------------------------------
 * Function:    matchMissmatchScore
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
int matchMissmatchScore(long long int i, long long int j, int8_t* a, int8_t* b) {
    if (a[i-1] == b[j-1])
        return matchScore;
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */