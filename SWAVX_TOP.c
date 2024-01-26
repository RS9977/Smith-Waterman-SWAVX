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
        m = 6000;
        n = 5000;
    }
    
    //Allocates a and b
    a = malloc((m+32) * sizeof(int8_t));
    b = malloc((n+32) * sizeof(int8_t));
    
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

    #if defined(SMP) && defined(L8)
    int8_t query_prof[32*m];
    
    for(int i=0; i<m; i++){
        memcpy(query_prof+i*32, iBlosum62_8bit + a[i]*32, 32 * sizeof(int8_t));
    }
    #else
    int8_t* query_prof;
    #endif

    #ifdef B512
    SWAVX_512_SeqToSeq_SubMat(a, b, H, P, m, n, NumOfTest, maxVal);
    #else
    SWAVX_256_SeqToSeq_SubMat(a, b, H, P, m, n, NumOfTest, maxVal, query_prof, 0);
    #endif
    
    //Gets final time
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double MeanTime = interval(time_start, time_stop)/NumOfTest;
    printf("Elapsed time: %f\n", MeanTime);
    printf("GCUPS: %f\n", (m-1)*(n-1)/(1e9*MeanTime));
    printf("maxVal: %d\n", *maxVal);

    //Frees similarity matrixes
    free(H);
    #ifdef BT
    free(P);
    #endif
    

    //Frees input arrays
    free(a);
    free(b);

    free(maxVal);
    
    return 0;

}
