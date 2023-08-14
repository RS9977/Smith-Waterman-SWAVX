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
    generate(a, b, m, n);

    SWAVX_256_SeqToSeq_SubMat(a, b, H, P, m, n, 1, -10);

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

}
