#include "SWAVX_TOP_datasets_MultiThread.h"

int main(int argc, char* argv[]) {
    
    double wk = wakeup_delay();

    const char *filenameA = "test2.fasta"; 
    const char *filenameB = "test3.fasta"; 
    int NumOfThreads      = getNumCPUThreads();
    if(NumOfThreads<1){
        NumOfThreads = 10;
    }
    if(argc==2)
    {
        NumOfThreads = strtoll(argv[1], NULL, 10);
    }
    else if(argc>2){
        filenameA    = argv[1];
        filenameB    = argv[2];
        NumOfThreads = strtoll(argv[3], NULL, 10);
    }
    

    //Reading datasets
    ProteinEntry *proteinEntriesA;
    int numEntriesA;
    if (!readProteinDataset(filenameA, &proteinEntriesA, &numEntriesA, 5e2)) {
        printf("Unable to read the file: %s\n",filenameA);
        return 0;        
    }
    

    ProteinEntry *proteinEntriesB;
    int numEntriesB;
    if (!readProteinDataset(filenameB, &proteinEntriesB, &numEntriesB, 1e7)) {
        printf("Unable to read the file: %s\n",filenameB);
        return 0;        
    }
    
    long long int HsizeA = 0;
    long long int HsizeB = 0;
    int i,j;
    
    int max_sizeA = 0;
    int max_SizeB = 0;

    for(i=0; i<numEntriesA; i++){
        HsizeA += proteinEntriesA[i].length;
        if (proteinEntriesA[i].length > max_sizeA)
            max_sizeA = proteinEntriesA[i].length;
    }

    for(i=0; i<numEntriesB; i++){
        HsizeB += proteinEntriesB[i].length;
        if (proteinEntriesB[i].length > max_SizeB)
            max_SizeB = proteinEntriesB[i].length;
    }
    
    //Load balancing
    int B_chunck_start[MAX_THREAD];
    int B_chunck_num  [MAX_THREAD];
    int B_chunck_size [MAX_THREAD];
    load_balance(B_chunck_start, B_chunck_num, B_chunck_size, HsizeB, numEntriesB, proteinEntriesB, NumOfThreads);

    int MaxHSize = (max_sizeA+2)*(max_SizeB+2);
    

    INT* maxVal  = calloc(numEntriesA * numEntriesB, sizeof(INT));
    //Allocates similarity matrix H
    //INT *H;
    #ifdef SAVEHP
    INT *H = calloc((HsizeA+numEntriesA) * (HsizeB+numEntriesB), sizeof(INT));
    #else
    INT *H = calloc(NumOfThreads*MaxHSize, sizeof(INT));
    #endif

    //Allocates predecessor matrix P
    
    #ifdef BT
    #ifdef SAVEHP
    INT *P = calloc((HsizeA+numEntriesA) * (HsizeB+numEntriesB), sizeof(INT));
    #else
    INT *P = calloc(NumOfThreads*MaxHSize, sizeof(INT));
    #endif
    #else
    INT *P = calloc(1, sizeof(INT));
    #endif

    struct timespec time_start, time_stop;
    clock_gettime(CLOCK_REALTIME, &time_start);
        for(int kk=0; kk<1; kk++){
        pthread_t threads[MAX_THREAD];
        WorkerIns thread_data_array[MAX_THREAD];
        int t;
        int rc;          
        long long int start = 0;
        for(t=0; t<NumOfThreads; t++){ 
            
            thread_data_array[t].proteinA = proteinEntriesA;  
            thread_data_array[t].proteinB = proteinEntriesB + B_chunck_start[t];  
            thread_data_array[t].H        = H+start;
            #ifdef BT
            thread_data_array[t].P        = P+start;
            #endif
            thread_data_array[t].A_num    = numEntriesA;
            thread_data_array[t].B_num    = B_chunck_num[t];
            thread_data_array[t].maxVal   = maxVal + numEntriesA* B_chunck_start[t];
            
            rc = pthread_create(&threads[t], NULL, chunck_computations, (void*) &thread_data_array[t]);
            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
            #ifdef SAVEHP
            start += (HsizeA+numEntriesA) * (B_chunck_size[t]+B_chunck_num[t]);
            #else
            start += MaxHSize;
            #endif
        }
        
        for (t = 0; t<NumOfThreads; t++) {
            if (pthread_join(threads[t], NULL)){ 
                printf("ERROR; code on return from join is %d\n", rc);
                exit(-1);
            }
        }
        }

    //Gets final time
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double MeanTime = interval(time_start, time_stop)/NumOfTest;
    printf("Elapsed time: %f\n", MeanTime);
    printf("GCUPS: %f\n", HsizeA*HsizeB/(1e9*MeanTime));
    printf("WakeUpVal: %f\n",wk);

    

    for(int i=0; i<numEntriesA; i++){
        printf("\nQuery number %d: \n", i);
        INT max[bestMatchNum]={0};
        int max_ind[bestMatchNum];
        for(int k=0; k<bestMatchNum; k++){
            for(int j=0; j<numEntriesB; j++){
                int ind = i+numEntriesA*j;
                if (!isValueInArray(max_ind, k, ind)){
                    if(max[k]<maxVal[ind]){
                        max[k]     = maxVal[ind];
                        max_ind[k] = ind;
                    }
                }
            }
        }
        for(int k=0; k<bestMatchNum; k++){
            printf("(j: %d, max: %d)| ",max_ind[k], max[k]);
        }
    }
    printf("\n");

    //Frees similarity matrixes
    free(H);
    free(P);

    free(maxVal);
    
    for (int i = 0; i < numEntriesA; i++) {
        free(proteinEntriesA[i].protein);
    }
    free(proteinEntriesA);

    for (int i = 0; i < numEntriesB; i++) {
        free(proteinEntriesB[i].protein);
    }
    free(proteinEntriesB);

    return 0;
}

//ٌWorker
void* chunck_computations(void* in){
    WorkerIns *inss = (WorkerIns *) in;
    ProteinEntry *proteinEntriesA = inss -> proteinA;
    ProteinEntry *proteinEntriesB = inss -> proteinB;
    INT* H                        = inss -> H;
    INT* maxVal                   = inss -> maxVal;
    INT* P;
    #ifdef BT
         P                        = inss -> P;
    #endif
    int  A_num                    = inss -> A_num;
    int  B_num                    = inss -> B_num;
    #ifdef SAVEHP
    long long int start           = 0;
    #endif
    int  i,j;
    
    for(i=0; i<A_num; i++){
        for(j=0; j<B_num; j++){
            #ifdef PARASAIL
                parasail_result_t *result = NULL;
                #ifdef P_scan
                    result = parasail_sw_scan_avx2_256_8(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #elif P_diag
                    result = parasail_sw_diag_avx2_256_8(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #else
                    result = parasail_sw_striped_avx2_256_32(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #endif
                *(maxVal+i+A_num*j) = result->score;
                parasail_result_free(result);
            #else
                #ifdef B512
                    #if SAVEHP
                        #ifdef BT
                        if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                            SWAVX_512_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j));
                        else
                            SWAVX_512_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j));
                        #else
                        if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                            SWAVX_512_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j));
                        else
                            SWAVX_512_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j));
                        #endif
                    #else
                    if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                        SWAVX_512_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j));
                    else
                        SWAVX_512_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j));
                    #endif
                #else
                    #ifdef SAVEHP
                        #ifdef BT
                        if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j));
                        else
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j));
                        #else
                        if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j));
                        else
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j));
                        #endif
                    #else
                    if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                        SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j));
                    else
                        SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j));
                    #endif
                #endif
            #endif
            #ifdef SAVEHP
            start += (proteinEntriesA[i].length+1) * (proteinEntriesB[j].length+1);
            #endif
        }
    }
    
}