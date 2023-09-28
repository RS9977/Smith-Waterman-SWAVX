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
    
    //Allocates a and b
    a = malloc(m * sizeof(int8_t)+8);
    b = malloc(n * sizeof(int8_t)+8);
    
    //Because now we have zeros
    m++;
    n++;
    
    
    
    
    INT *H;
    H = malloc(m * n* sizeof(INT));

    *(H) = 0;
    *(H+1) = 0;
    *(H+2) = 0;


    
    //Allocates predecessor matrix P
    INT *P;
    #ifdef BT
    P = calloc(m * n, sizeof(INT));
    #endif

    
    

    //Gen rand arrays a and b
    generate(a, b, m, n);

    struct timespec time_start, time_stop;
    clock_gettime(CLOCK_REALTIME, &time_start);
    
    INT* maxVal = calloc(1, sizeof(INT));

    #ifdef B512
    SWAVX_512_SeqToSeq_SubMat(a, b, H, P, m, n, NumOfTest, maxVal);
    #else
    SWAVX_256_SeqToSeq_SubMat(a, b, H, P, m, n, NumOfTest, maxVal);
    #endif
    
    //Gets final time
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double MeanTime = interval(time_start, time_stop)/NumOfTest;
    printf("Elapsed time: %f\n", MeanTime);
    printf("GCUPS: %f\n", (m-1)*(n-1)/(1e9*MeanTime));
    printf("maxVal: %d\n", *maxVal);

    //Frees similarity matrixes
    free(H);
    free(P);

    //Frees input arrays
    free(a);
    free(b);

    free(maxVal);

    return 0;

}
