/***********************************************************************
 * Smith–Waterman algorithm
 * Purpose:     Local alignment of nucleotide or protein sequences
 * Authors:     Daniel Holanda, Hanoch Griner, Taynara Pinheiro
 ***********************************************************************/
//This is the working AVX256 for 4B integer

//gcc -mavx2 -O3 SWAVX_256_SingleThread_SubMat.c -lgomp -o SWAVX_256_SingleThread_SubMat
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

//#include <zmmintrin.h>


int const iBlosum62 [] = {
5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5, -3, 
-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 5 
};

/*--------------------------------------------------------------------
 * Text Tweaks
 */
#define RESET   "\033[0m"
#define BOLDRED "\033[1m\033[31m"      /* Bold Red */
/* End of text tweaks */

/*--------------------------------------------------------------------
 * Constants
 */
#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

#define Vsize 8


#define NumOfTest 1e0//1e4
//#define DEBUG
//#define pragmas
/* End of constants */


/*--------------------------------------------------------------------
 * Functions Prototypes
 */
void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, int* H, int* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len);
void similarityScoreIntrinsic(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, long long int ii, long long int jj, int* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len);
int  matchMissmatchScore(long long int i, long long int j);
void backtrack(int* P, long long int maxPos, long long int maxPos_max_len);
void printMatrix(int* matrix);
void printPredecessorMatrix(int* matrix);
void generate(void);
/* End of prototypes */


/*--------------------------------------------------------------------
 * Global Variables
 */
//Defines size of strings to be compared
long long int m = 11; //Columns - Size of string a
long long int n = 7;  //Lines - Size of string b

//Defines scores
int matchScore = 5;
int missmatchScore = -3;
int gapScore = -10; 

//Strings over the Alphabet Sigma
int8_t *a, *b;

/* End of global variables */

/*--------------------------------------------------------------------
 * Function:    main
 */
int main(int argc, char* argv[]) {
    
    
    if(argc>1){
        m = strtoll(argv[1], NULL, 10);
        n = strtoll(argv[2], NULL, 10); 
        long long int temp;
        if( m<n){
            temp = m;
            m = n;
            n = temp;
        }
    }
    else{
        m = 9;
        n = 8;
    }


    #ifdef DEBUG
    printf("\nMatrix[%lld][%lld]\n", n, m);
    #endif


    
    //Allocates a and b
    a = malloc(m * sizeof(int8_t));
    b = malloc(n * sizeof(int8_t));
    
    //Because now we have zeros
    m++;
    n++;
    
    //Allocates similarity matrix H
    int *H;
    H = calloc(m * n, sizeof(int));

    //Allocates predecessor matrix P
    int *P;
    P = calloc(m * n, sizeof(int));


    //Gen rand arrays a and b
    generate();
    // a[0] =   'C';
    // a[1] =   'G';
    // a[2] =   'T';
    // a[3] =   'G';
    // a[4] =   'A';
    // a[5] =   'A';
    // a[6] =   'T';
    // a[7] =   'T';
    // a[8] =   'C';
    // a[9] =   'A';
    // a[10] =  'T';

    // b[0] =   'G';
    // b[1] =   'A';
    // b[2] =   'C';
    // b[3] =   'T';
    // b[4] =   'T';
    // b[5] =   'A';
    // b[6] =   'C';



    //Start position for backtrack
    long long int maxPos         = 0;
    long long int maxPos_max_len = 0;

    //Calculates the similarity matrix
    long long int i, j;

    
    double initialTime = omp_get_wtime();
    #ifdef pragmas
    #pragma GCC ivdep
    #endif


    #ifdef DEBUG
    printf("\n a string:\n");
    for(i=0; i<m-1; i++)
        printf("%d ",a[i]);
    printf("\n b string:\n");
    for(i=0; i<n-1; i++)
        printf("%d ",b[i]);
    printf("\n");
    #endif
    double t;
    int it;
    for(it=0; it<NumOfTest; it++){

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
        }
	    else if (i>=m){
		    max_len = m+n-1-i;
            j_start = 0;   
            j_end   = max_len;     
            ind_u   = ind - max_len - 1;
            ind_l   = ind - max_len; 
            ind_d   = ind - (max_len<<1) - 2;
        }
	    else{
		    max_len   = n;
            j_start = 1;
            j_end   = max_len;
            ind_u     = ind - max_len - 1;
            ind_l     = ind - max_len; 
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
     
        for (j = j_start; j <j_end-Vsize+1; j+=Vsize) { //Columns  

           if (i<m){
                ii = i-j;
                jj = j;
            }
            else{
                ii = m-1-j;
                jj = i-m+j+1;
            } 
           similarityScoreIntrinsic(HH, Hu, Hd, Hl, PP, reverseIndices, ii, jj, H, ind+j, max_len, &maxPos, &maxPos_max_len);
           Hu++;
           Hl++;
           Hd++;
           HH++;
           PP++;
        }
       
        for(;j<j_end; j++){
            if (i<m){
                ii = i-j;
                jj = j;
            }
            else{
                ii = m-1-j;
                jj = i-m+j+1;
            }      
            similarityScore(ind+j, ind_u+j, ind_d+j, ind_l+j, ii, jj, H, P, max_len, &maxPos, &maxPos_max_len);
        }
        ind += max_len;
    }
    
    backtrack(P, maxPos, maxPos_max_len);


    }


    //Gets final time
    double finalTime = omp_get_wtime();

    printf("\nElapsed time: %f\n\n", (finalTime - initialTime)/NumOfTest);

    #ifdef DEBUG
    FILE* file = fopen("256_ST_SB.txt", "w");

    // Check if the file was opened successfully
    if (file == NULL) {
        printf("Failed to open the file.\n");
        return 1; // Return an error code
    }
    int ind;
    fprintf(file, " \t \t");
    for(i=0; i<m-1; i++)
        fprintf(file,"%d\t",a[i]);
    fprintf(file, "\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                fprintf(file, " \t");
            else if(j==-1 && i>0)
                fprintf(file, "%d\t",b[i-1]); 
            else
                fprintf(file, "%d\t", H[ind]);
        }
        fprintf(file, "\n");
    }
    #endif
    
    #ifdef DEBUG
    printf("\nSimilarity Matrix:\n");
    printMatrix(H);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P);
    #endif

    //Frees similarity matrixes
    free(H);
    free(P);

    //Frees input arrays
    free(a);
    free(b);

    return 0;
}  /* End of main */


/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate  the maximum Similarity-Score H(i,j)
 */
void similarityScore(long long int ind, long long int ind_u, long long int ind_d, long long int ind_l, long long int ii, long long int jj, int* H, int* P, long long int max_len, long long int* maxPos, long long int *maxPos_max_len) {

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


void similarityScoreIntrinsic(__m256i* HH,__m256i* Hu,__m256i* Hd,__m256i* Hl,__m256i* PP, __m256i reverseIndices, long long int ii, long long int jj, int* H, long long int ind, long long int max_len, long long int* maxPos, long long int *maxPos_max_len) {

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
 * Function:    matchMissmatchScore
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
int matchMissmatchScore(long long int i, long long int j) {
    if (a[i-1] == b[j-1])
        return matchScore;
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */

/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(int* P, long long int maxPos, long long int maxPos_max_len) {
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
 * Function:    printMatrix
 * Purpose:     Print Matrix
 */
void printMatrix(int* matrix) {
    long long int i, j, ind;
    printf(" \t \t");
    for(i=0; i<m-1; i++)
        printf("%d\t",a[i]);
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                printf(" \t");
            else if(j==-1 && i>0)
                printf("%d\t",b[i-1]); 
            else
                printf("%d\t", matrix[ind]);
        }
        printf("\n");
    }
}  /* End of printMatrix */

/*--------------------------------------------------------------------
 * Function:    printPredecessorMatrix
 * Purpose:     Print predecessor matrix
 */
void printPredecessorMatrix(int* matrix) {
    long long int i, j, ind;
    printf("    ");
    for(i=0; i<m-1; i++)
        printf("%d ",a[i]);
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                printf("  ");
            else if(j==-1 && i>0)
                printf("%d ",b[i-1]); 
            else{
                if(matrix[ind] < 0) {
                    printf(BOLDRED);
                    if (matrix[ind] == -UP)
                        printf("↑ ");
                    else if (matrix[ind] == -LEFT)
                        printf("← ");
                    else if (matrix[ind] == -DIAGONAL)
                        printf("↖ ");
                    else
                        printf("- ");
                    printf(RESET);
                } else {
                    if (matrix[ind] == UP)
                        printf("↑ ");
                    else if (matrix[ind] == LEFT)
                        printf("← ");
                    else if (matrix[ind] == DIAGONAL)
                        printf("↖ ");
                    else
                        printf("- ");
            }
            }
        }
        printf("\n");
    }

}  /* End of printPredecessorMatrix */


/*--------------------------------------------------------------------
 * Function:    generate
 * Purpose:     Generate arrays a and b
 */
 void generate(){
    //Generates the values of a
    long long int i;
    for(i=0;i<m;i++){
        a[i] = rand()%4;
        /*
        int aux=rand()%4;
        if(aux==0)
            a[i]='A';
        else if(aux==2)
            a[i]='C';
        else if(aux==3)
            a[i]='G';
        else
            a[i]='T';
        */
    }

    //Generates the values of b
    for(i=0;i<n;i++){
        b[i] = rand()%4;
        /*
        int aux=rand()%4;
        if(aux==0)
            b[i]='A';
        else if(aux==2)
            b[i]='C';
        else if(aux==3)
            b[i]='G';
        else
            b[i]='T';
        */
    }
} /* End of generate */


/*--------------------------------------------------------------------
 * External References:
 * http://vlab.amrita.edu/?sub=3&brch=274&sim=1433&cnt=1
 * http://pt.slideshare.net/avrilcoghlan/the-smith-waterman-algorithm
 * http://baba.sourceforge.net/
 */