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

    double initialTime = omp_get_wtime();

    #ifdef DEBUG
    printf("\n a string:\n");
    for(i=0; i<m-1; i++)
        printf("%d ",a[i]);
    printf("\n b string:\n");
    for(i=0; i<n-1; i++)
        printf("%d ",b[i]);
    printf("\n");
    #endif
    
    SWAVX_256_SeqToSeq_SubMat(a, b, H, P, m, n, NumOfTest, -10);
    
    //Gets final time
    double finalTime = omp_get_wtime();

    double MeanTime = (finalTime - initialTime)/NumOfTest;
    printf("\nElapsed time: %f\n", MeanTime);
    printf("GCUPS: %f\n", (m-1)*(n-1)/(1e9*MeanTime));


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
