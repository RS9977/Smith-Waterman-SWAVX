#include "SWAVX_TOP.h"

int main(int argc, char* argv[]) {
    
    int m, n;
    int8_t *a, *b;
    
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
    printf("\nMatrix[%d][%d]\n", n, m);
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
    generate(a, b, m, n);

    SWAVX_256_SeqToSeq_SubMat(a, b, H, P, m, n, 1, -10);
    
    #ifdef DEBUG
    saveInFile(H, a, b, m, n);

    printf("\nSimilarity Matrix:\n");
    printMatrix(H, a, b, m, n);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P, a, b, m, n);
    #endif

    //Frees similarity matrixes
    free(H);
    free(P);

    //Frees input arrays
    free(a);
    free(b);

    return 0;

}
