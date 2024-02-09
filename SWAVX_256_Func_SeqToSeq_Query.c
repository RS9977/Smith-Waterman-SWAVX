#include "SWAVX_256_Func_SeqToSeq_Query.h"
_Thread_local __m256i local_max;
_Thread_local __m256i scores_mat[32];
void SWAVX_256_SeqToSeq_QueryLTBatch(ProteinBatch batch, int8_t* b, INT *H, INT* P, int n, int NumOfTest, INT* maxVal, int8_t* query_prof, int MaxHSize){
    
    if(MaxHSize>BLOCKSIZE){
    int blockNum = MaxHSize/BLOCKSIZE;
    MaxHSize /= blockNum;
    }
    

    int m = batch.lengths[31];
    //Calculates the similarity matrix
    int i, j;

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
    
    
    __m256i* query_prof_pointer = (__m256i*) query_prof;
    local_max = _mm256_setzero_si256();    

    //__m256i *local_max = malloc(sizeof(__m256i));
    
    #ifdef L8
   
    int Vsize = 32;

    
    #elif L16
  

    int Vsize = 16;

    #else
    //int* local_max_mem = malloc(8*sizeof(int));
    //__m256i *local_max = (__m256i*) local_max_mem;
    
  
    
    int Vsize = 8;
    #endif
   // _mm256_storeu_si256(local_max, temp);
    
    int score_cnt = 0;
    
    double t;
    int it;
    
    int maxPos         = 0;
    int maxPos_max_len = 0;
    

    int ind   = 3;
    int indd  = 0;
    int indul = 1;
    int ind_u, ind_d, ind_l; 
    int block_cnt = 0;
   __m256i offset =_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    
    
    for (i = 2; i < m+n-1; i++) { //Lines
        int max_len;
        int ii,jj;
        int j_start, j_end;
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
        
        #ifndef AFFINE
        /*
        __m256i* Hu = (__m256i*) (H+ind_u+j_start+k*MaxHSize);
        __m256i* Hl = (__m256i*) (H+ind_l+j_start+k*MaxHSize);
        __m256i* Hd = (__m256i*) (H+ind_d+j_start+k*MaxHSize);
        __m256i* HH = (__m256i*) (H+ind+j_start+k*MaxHSize);
        __m256i* PP; 
        */
        #ifdef BT
        PP = (__m256i*) (P+ind+j_start+k*MaxHSize);
        #endif
        for (j=j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns 
            int score_cnt = 0;
            #pragma unroll 32
            for(int k=0; k<32; k++){
                int bias = j_start+k*MaxHSize+(j-j_start)+block_cnt*BLOCKSIZE*32;
                __m256i* Hu = (__m256i*) (H+ind_u+bias);
                __m256i* Hl = (__m256i*) (H+ind_l+bias);
                __m256i* Hd = (__m256i*) (H+ind_d+bias);
                __m256i* HH = (__m256i*) (H+ind+bias);
                __m256i* PP;  
            if(k==0){ 
            #pragma unroll 32
                for(int p=0; p<32; p++){
                    __m256i* query_prof_pointer = (__m256i*) query_prof;
                    __m256i query_temp = _mm256_loadu_si256(query_prof_pointer+p+jj-1);
                    __m256i b_temp = _mm256_loadu_si256((__m256i*)(batch.transposedSequences+p+jj-1));
                    __m256i temp = _mm256_shuffle_epi8(query_temp, b_temp);
                    _mm256_storeu_si256(scores_mat+p+score_cnt*32, temp);
                }
                
                //#pragma GCC inline
                Transpose32By32(scores_mat+score_cnt*32);
                score_cnt++;
            }
            
                
                #pragma GCC inline      
                similarityScoreIntrinsic8(HH, Hu, Hd, Hl, PP, ii, jj, H, ind+j, max_len, &maxPos, &maxPos_max_len, maxVal, k, (score_cnt-1)*32, 32);
                }
                ii -= Vsize;
                jj += Vsize;
                /*
                Hu++;
                Hl++;
                Hd++;
                HH++;
                PP++;
                */
            }
            
            int score_cnt = 0;
            #pragma unroll 32
            for(int k=0; k<32; k++){
                int bias = j_start+k*MaxHSize+(j-j_start)+block_cnt*BLOCKSIZE*32;
                __m256i* Hu = (__m256i*) (H+ind_u+bias);
                __m256i* Hl = (__m256i*) (H+ind_l+bias);
                __m256i* Hd = (__m256i*) (H+ind_d+bias);
                __m256i* HH = (__m256i*) (H+ind+bias);
                __m256i* PP; 
            
            if(j_end-j > 4){
                if(k==0){ 
            #pragma unroll 32
                for(int p=0; p<32; p++){
                    __m256i* query_prof_pointer = (__m256i*) query_prof;
                    __m256i query_temp = _mm256_loadu_si256(query_prof_pointer+p+jj-1);
                    __m256i b_temp = _mm256_loadu_si256((__m256i*)(batch.transposedSequences+p+jj-1));
                    __m256i temp = _mm256_shuffle_epi8(query_temp, b_temp);
                    _mm256_storeu_si256(scores_mat+p+score_cnt*32, temp);
                }
                //#pragma GCC inline
                
                Transpose32By32(scores_mat+score_cnt*32);
                score_cnt++;
            }

                similarityScoreIntrinsic8(HH, Hu, Hd, Hl, PP, ii, jj, H, ind+j, max_len, &maxPos, &maxPos_max_len, maxVal, k, (score_cnt-1)*32, j_end-j-1);               
                
            }
            }
            /*
            for(int k=0; k<32; k++){
            int j_temp = jj;
            int i_temp = ii;
            for(;j<j_end; j++){
                #pragma GCC inline   
               // similarityScore_QueryLTBatch(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H+k*MaxHSize, P, max_len, &maxPos, &maxPos_max_len, maxVal, batch, b, n, k);
                ii --;
                jj ++;
            }
            jj = j_temp;
            ii = i_temp;
            }
            */
            ind += max_len;
            if(ind>BLOCKSIZE*(block_cnt+1)){
                block_cnt ++;
            }

        }
        
        #else
        
        
        #endif
   
    
    
       #ifdef BT
       backtrack(P, maxPos, maxPos_max_len, m, n);
       #endif
    

    #ifndef BT
    __m256i max = local_max;
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

   // free(local_max);
    
}
void SWAVX_256_SeqToSeq_QueryBTBatch(int8_t *a, ProteinBatch batch, INT *H, INT* P, int m, int NumOfTest, INT* maxVal, int8_t* query_prof, int MaxHSize){


    
    //Calculates the similarity matrix
    int i, j;
    int n = batch.lengths[31];
    
    if(MaxHSize>BLOCKSIZE){
    int blockNum = MaxHSize/BLOCKSIZE;
    MaxHSize /= blockNum;
    }

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
    
    
    __m256i* query_prof_pointer = (__m256i*) query_prof;
    local_max = _mm256_setzero_si256();


    
    
    

    //__m256i *local_max = malloc(sizeof(__m256i));
    
    #ifdef L8
   
    int Vsize = 32;
    #ifdef SUBMAT

    #else
   
    #endif
    
    #elif L16
  

    int Vsize = 16;
    #else
    //int* local_max_mem = malloc(8*sizeof(int));
    //__m256i *local_max = (__m256i*) local_max_mem;
    
  
    
    int Vsize = 8;

    #endif
   // _mm256_storeu_si256(local_max, temp);
   
    
    double t;
    int it;
    
    int maxPos         = 0;
    int maxPos_max_len = 0;
    

    int ind   = 3;
    int indd  = 0;
    int indul = 1;
    int ind_u, ind_d, ind_l; 
   __m256i offset =_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    
    int block_cnt = 0;
    for (i = 2; i < m+n-1; i++) { //Lines
        int max_len;
        int ii,jj;
        int j_start, j_end;
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

        #ifndef AFFINE
       
        for (j=j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns  
            int score_cnt = 0;
                #pragma unroll 32
                for(int k=0; k<32; k++){
                int bias = j_start+k*MaxHSize+(j-j_start) + block_cnt*BLOCKSIZE*32;
                __m256i* Hu = (__m256i*) (H+ind_u+bias);
                __m256i* Hl = (__m256i*) (H+ind_l+bias);
                __m256i* Hd = (__m256i*) (H+ind_d+bias);
                __m256i* HH = (__m256i*) (H+ind+bias);
                __m256i* PP;  
                if(k==0){
                    #pragma unroll 32
                    for(int p=0; p<32; p++){
                        __m256i* query_prof_pointer = (__m256i*) query_prof;
                        __m256i query_temp = _mm256_loadu_si256(query_prof_pointer+p+ii-1);
                        __m256i b_temp = _mm256_loadu_si256((__m256i*)(batch.transposedSequences+p+ii-1));               
                        __m256i temp = _mm256_shuffle_epi8(query_temp, b_temp);
                        _mm256_storeu_si256(scores_mat+p+score_cnt*32, temp);
                        
                    }
                    
                    #pragma GCC inline
                    Transpose32By32(scores_mat+score_cnt*32);
                    score_cnt++;
                }
                #ifdef BT
                PP = (__m256i*) (P+ind+j_start+k*MaxHSize);
                #endif
                #pragma GCC inline      
                similarityScoreIntrinsic8(HH, Hu, Hd, Hl, PP, ii, jj, H, ind+j, max_len, &maxPos, &maxPos_max_len, maxVal, k, (score_cnt-1)*32, 32);               
                }
                ii -= Vsize;
                jj += Vsize;
                /*
                Hu++;
                Hl++;
                Hd++;
                HH++;
                PP++;
                */
            }
            
            int score_cnt = 0;
            #pragma unroll 32
            for(int k=0; k<32; k++){
            int bias = j_start+k*MaxHSize+(j-j_start)+block_cnt*BLOCKSIZE*32;
                __m256i* Hu = (__m256i*) (H+ind_u+bias);
                __m256i* Hl = (__m256i*) (H+ind_l+bias);
                __m256i* Hd = (__m256i*) (H+ind_d+bias);
                __m256i* HH = (__m256i*) (H+ind+bias);
                __m256i* PP; 
            if(j_end-j > 4){
                if(k==0){
                    #pragma unroll 32
                    for(int p=0; p<32; p++){
                        __m256i* query_prof_pointer = (__m256i*) query_prof;
                        __m256i query_temp = _mm256_loadu_si256(query_prof_pointer+p+ii-1);
                        __m256i b_temp = _mm256_loadu_si256((__m256i*)(batch.transposedSequences+p+ii-1));               
                        __m256i temp = _mm256_shuffle_epi8(query_temp, b_temp);
                        _mm256_storeu_si256(scores_mat+p+score_cnt*32, temp);
                    }
                    
                    #pragma GCC inline
                    Transpose32By32(scores_mat+score_cnt*32);
                    score_cnt++;
                }
                
                #pragma GCC inline
                similarityScoreIntrinsic8(HH, Hu, Hd, Hl, PP, ii, jj, H, ind+j, max_len, &maxPos, &maxPos_max_len, maxVal, k, (score_cnt-1)*32, j_end-j-1);               
                j = j_end;
            }
            }
            /*
            for(int k=0; k<32; k++){
            int ii_temp = ii;
            int jj_temp = jj;
            for(;j<j_end; j++){   
                #pragma GCC inline
             //   similarityScore_QueryBTBatch(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H+k*MaxHSize, P, max_len, &maxPos, &maxPos_max_len, maxVal, a, batch, m, k);
                ii --;
                jj ++;
            }
            ii = ii_temp;
            jj = jj_temp;
            }
            */
            ind += max_len;
            if(ind>BLOCKSIZE*(block_cnt+1)){
                block_cnt ++;
               // printf("block size %d\n", block_cnt);
            }
            
        }       
        #else
        
        
        #endif 
    
       #ifdef BT
       backtrack(P, maxPos, maxPos_max_len, m, n);
       #endif
    

    #ifndef BT
    __m256i max = local_max;
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

   // free(local_max);
   
   
}

void similarityScore_QueryLTBatch(int ind, int ind_u, int ind_d, int ind_l, int ii, int jj, INT* H, INT* P, int max_len, int* maxPos, int *maxPos_max_len, INT* maxVal, ProteinBatch batch, int8_t *b, int n, int k) {

    INT up, left, diag;

    //Get element above
    up   = H[ind_u] + gapScore;

    //Get element on the left
    left = H[ind_l] + gapScore;


    //Get element on the diagonal
    diag = H[ind_d] + iBlosum62[batch.transposedSequences[(ii-1)*32+k]*32+b[jj-1]];


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

void similarityScore_QueryBTBatch(int ind, int ind_u, int ind_d, int ind_l, int ii, int jj, INT* H, INT* P, int max_len, int* maxPos, int *maxPos_max_len, INT* maxVal, int8_t *a, ProteinBatch batch, int m, int k) {

    INT up, left, diag;

    //Get element above
    up   = H[ind_u] + gapScore;

    //Get element on the left
    left = H[ind_l] + gapScore;

  
    diag = H[ind_d] + iBlosum62[a[ii-1]*32+batch.transposedSequences[(jj-1)*32+k]];
    

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



void similarityScoreIntrinsic8(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, int ii, int jj, INT* H, int ind, int max_len, int* maxPos, int *maxPos_max_len, INT* maxVal, int k, int score_bias, int reminder) {
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

    diag = _mm256_add_epi8(HHd, _mm256_loadu_si256(scores_mat+k+score_bias));

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
    //__m256i tempmax = _mm256_loadu_si256(local_max);
    local_max = _mm256_max_epi8(max, local_max);
    //_mm256_storeu_si256(local_max,tempmax);
    #endif
}  /* End of similarityScore */



/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(INT* P, int maxPos, int maxPos_max_len, int m, int n) {
    //hold maxPos value
    int predPos;
    int predMaxLen;
    #ifdef pragmas
    #pragma GCC ivdep
    #endif
    //backtrack from maxPos to startPos = 0 
    int first_sec = (n*(n+1))/2;
    int last_sec  = n*m - (n*(n-1))/2;
    int ind_u, ind_d, ind_l; 
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
int matchMissmatchScore(int i, int j, int8_t* a, int8_t* b) {
    if (a[i-1] == b[j-1])
        return matchScore;
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */

void Transpose32By32(__m256i mat[32]) {
    __m256i decode_mask = _mm256_setr_epi8(0xF0, 0x00, 0x00, 0x00,
                                               0x0F, 0x00, 0x00, 0x00,
                                               0x00, 0xF0, 0x00, 0x00,
                                               0x00, 0x0F, 0x00, 0x00,
                                               0x00, 0x00, 0xF0, 0x00,
                                               0x00, 0x00, 0x0F, 0x00,
                                               0x00, 0x00, 0x00, 0xF0,
                                               0x00, 0x00, 0x00, 0x0F);
    __m256i decode_sub = _mm256_setr_epi32(1<<28-1, 1<<24-1, 1<<20-1, 1<<16-1, 1<<12-1, 1<<8-1, 1<<4-1, 0);
    __m256i temp2;
    __m256i temp; // Temporary registers for holding the data
    __m256i accumulators[4] = {
        _mm256_setzero_si256(),
        _mm256_setzero_si256(),
        _mm256_setzero_si256(),
        _mm256_setzero_si256()
    }; // Initialize accumulators with zero

    __m256i transposed_32[128];
    #pragma unroll 4
    for (int i = 0; i < 4; i ++) { // Process 8 rows at a time
        #pragma unroll 8
        for (int j = 0; j < 8; ++j) {
            temp = _mm256_loadu_si256((__m256i*)&mat[i*8+j]); 
            __m256i zero = _mm256_setzero_si256();

            // Split the 8-bit integers into two groups of 16, each expanded to 16 bits
            __m256i lo16 = _mm256_unpacklo_epi8(temp, zero); // Low 16 bytes of the row, expanded
            __m256i hi16 = _mm256_unpackhi_epi8(temp, zero); // High 16 bytes of the row, expanded

            // Now, split each of those into 32-bit integers
            __m256i row0 = _mm256_unpackhi_epi16(hi16, zero); 
            __m256i row1 = _mm256_unpacklo_epi16(lo16, zero); // Low 16 bytes -> low 8 integers
            __m256i row2 = _mm256_unpackhi_epi16(lo16, zero); // Low 16 bytes -> high 8 integers
            __m256i row3 = _mm256_unpacklo_epi16(hi16, zero); // High 16 bytes -> low 8 integers


            if(j!=0){
                accumulators[0] = _mm256_slli_epi32(row0,4);
                accumulators[1] = _mm256_slli_epi32(row1,4);
                accumulators[2] = _mm256_slli_epi32(row2,4);
                accumulators[3] = _mm256_slli_epi32(row3,4);
            }
            accumulators[0] = _mm256_add_epi32(row0, accumulators[0]);
            accumulators[1] = _mm256_add_epi32(row1, accumulators[1]);
            accumulators[2] = _mm256_add_epi32(row2, accumulators[2]);
            accumulators[3] = _mm256_add_epi32(row3, accumulators[3]);
        }
        
        
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],0));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+0] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],1));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+1] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],2));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+2] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],3));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+3] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],4));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+4] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],5));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+5] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],6));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+6] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[0],7));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+7] = _mm256_sub_epi32(temp2, decode_sub);


        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],0));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+8] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],1));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+9] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],2));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+10] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],3));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+11] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],4));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+12] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],5));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+13] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],6));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+14] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[1],7));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+15] = _mm256_sub_epi32(temp2, decode_sub);



        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],0));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+16] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],1));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+17] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],2));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+18] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],3));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+19] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],4));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+20] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],5));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+21] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],6));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+22] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[2],7));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+23] = _mm256_sub_epi32(temp2, decode_sub);


        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],0));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+24] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],1));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+25] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],2));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+26] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],3));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+27] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],4));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+28] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],5));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+29] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],6));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+30] = _mm256_sub_epi32(temp2, decode_sub);
        temp2 = _mm256_set1_epi32(_mm256_extract_epi32(accumulators[3],7));
        temp2 = _mm256_and_si256(temp2, decode_mask);
        transposed_32[i*8+31] = _mm256_sub_epi32(temp2, decode_sub);

    }
    #pragma unroll 4
    for (int i = 0; i < 4; i ++) { // Process 8 rows at a time
        #pragma unroll 8
        for (int j = 0; j < 8; ++j) {
        __m256i gatheredData1     = _mm256_packs_epi32(transposed_32[i*8+j+0], transposed_32[i*8+j+8]);
        gatheredData1             = _mm256_permute4x64_epi64(gatheredData1, 0xd8);
        __m256i gatheredData2     = _mm256_packs_epi32(transposed_32[i*8+j+16], transposed_32[i*8+j+24]);
        gatheredData2             = _mm256_permute4x64_epi64(gatheredData2, 0xd8); 
        __m256i gatheredData      = _mm256_packs_epi16(gatheredData1,gatheredData2);
        gatheredData             = _mm256_permute4x64_epi64(gatheredData, 0xd8);
        _mm256_storeu_si256((__m256i*)&mat[i*8+j], gatheredData);
        }
    }
}