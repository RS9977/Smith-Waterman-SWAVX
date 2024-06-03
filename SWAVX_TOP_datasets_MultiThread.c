#include "SWAVX_TOP_datasets_MultiThread.h"

int main(int argc, char* argv[]) {
    
    //double wk = wakeup_delay();

    const char *filenameA = "test14.fasta"; 
    const char *filenameB = "test7.fasta"; 
    //int NumOfThreads      = getNumCPUThreads();
    int NumOfThreads      = 1;
    if(NumOfThreads<1){
        NumOfThreads = 1;
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
    
    printf("size B: %lld\n",  HsizeB);
    //Load balancing
    int B_chunck_start[MAX_THREAD];
    int B_chunck_num  [MAX_THREAD];
    int B_chunck_size [MAX_THREAD];
    #ifndef Query
    load_balance(B_chunck_start, B_chunck_num, B_chunck_size, HsizeB, numEntriesB, proteinEntriesB, NumOfThreads);
    #endif

    #ifndef Query
    int MaxHSize = (max_sizeA+2)*(max_SizeB+2);
    #else
    int MaxHSize = 32*(max_sizeA+2)*(max_SizeB+2);
    #endif

    ProteinBatch* batches;
    #ifdef Query
    int batchSizes = 32;
    int totalBatches = (numEntriesB + batchSizes - 1) / batchSizes; // Calculate the total number of batches, rounding up

    // Dynamically allocate memory for the batches array based on the calculated total number of batches
    batches = (ProteinBatch*)malloc(totalBatches * sizeof(ProteinBatch));
    if (batches == NULL) {
        fprintf(stderr, "Memory allocation for batches failed.\n");
        return 1;
    }
    for (int i = 0; i < totalBatches; ++i) {
        batches[i] = transposeProteins(&proteinEntriesB[i * batchSizes], batchSizes);
    }
    
    load_balance_batch(B_chunck_start, B_chunck_num, B_chunck_size, HsizeB, totalBatches, batches, NumOfThreads, 32);
    
    #endif

    INT* maxVal  = calloc(numEntriesA * numEntriesB, sizeof(INT));
    //Allocates similarity matrix H
    //INT *H;
    #ifndef DAlloc
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
    #else
    INT *H;
    #endif

    struct timespec time_start, time_stop;
    clock_gettime(CLOCK_REALTIME, &time_start);
        for(int kk=0; kk<NumOfTest; kk++){
        pthread_t threads[MAX_THREAD];
        WorkerIns thread_data_array[MAX_THREAD];
        int t;
        int rc;          
        long long int start = 0;
        for(t=0; t<NumOfThreads; t++){ 
            
            thread_data_array[t].proteinA = proteinEntriesA;  
            #ifndef Query
            thread_data_array[t].proteinB = proteinEntriesB + B_chunck_start[t];  
            
            #else
            thread_data_array[t].batch    = batches + B_chunck_start[t];
            thread_data_array[t].MaxHSize = MaxHSize;
            #endif
            #ifndef DAlloc
            thread_data_array[t].H        = H+start;
            #ifdef BT
            thread_data_array[t].P        = P+start;
            #endif
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
    //printf("WakeUpVal: %f\n",wk);

    

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
    #ifdef BT
    free(P);
    #endif

    free(maxVal);
    
    for (int i = 0; i < numEntriesA; i++) {
        free(proteinEntriesA[i].protein);
    }
    free(proteinEntriesA);

    for (int i = 0; i < numEntriesB; i++) {
        free(proteinEntriesB[i].protein);
    }
    free(proteinEntriesB);

    #ifdef Query
    for (int i = 0; i < totalBatches; ++i) {
        free(batches[i].transposedSequences);
        free(batches[i].lengths);
    }
    free(batches);
    #endif
    return 0;
}

//ٌWorker
void* chunck_computations(void* in){
    WorkerIns *inss = (WorkerIns *) in;
    ProteinEntry *proteinEntriesA = inss -> proteinA;
    #ifndef Query
    ProteinEntry *proteinEntriesB = inss -> proteinB;
   
    #else
    ProteinBatch *batches         = inss -> batch;
    int MaxHSize                  = inss -> MaxHSize;
    #endif
    #ifndef DAlloc
    INT* H                        = inss -> H;
    #endif
    
    INT* P;
    
    INT* maxVal                   = inss -> maxVal;
    #ifdef BT
         P                        = inss -> P;
    #endif
    int  A_num                    = inss -> A_num;
    int  B_num                    = inss -> B_num;
    #ifdef SAVEHP
    long long int start           = 0;
    #endif
    

    #ifdef PARASAIL
    parasail_matrix_t *matrix = parasail_matrix_create("ARNDCEQGHILKMFPSTWYV", matchScore, missmatchScore);
    #endif

    int  i,j;
    
    
    
    for(i=0; i<A_num; i++){
        #if defined(SMP) && defined(L8)
        int8_t query_prof[32*proteinEntriesA[i].length];    
        for(int i=0; i<proteinEntriesA[i].length; i++){
            memcpy(query_prof+i*32, iBlosum62_8bit + proteinEntriesA[i].protein[i]*32, 32 * sizeof(int8_t));
        }
        #elif defined(SMPP) && defined(L8)
        int8_t query_prof[32*proteinEntriesA[i].length];    
        for(int i=0; i<proteinEntriesA[i].length; i++){
            memcpy(query_prof+i*32, iBlosum62_8bit + proteinEntriesA[i].protein[i]*32, 32 * sizeof(int8_t));
        }
        #elif Query
        int8_t query_prof[32*proteinEntriesA[i].length];    
        for(int i=0; i<proteinEntriesA[i].length; i++){
            memcpy(query_prof+i*32, iBlosum62_8bit + proteinEntriesA[i].protein[i]*32, 32 * sizeof(int8_t));
        }
        #else
        int8_t* query_prof;
        #endif
        
        #ifdef PARASAIL
        #ifdef QP
        // Create a profile for the query sequence
        parasail_profile_t *profile = parasail_profile_create_avx_256_8(proteinEntriesA[i].protein, proteinEntriesA[i].length, &parasail_blosum62);
        #endif
        #endif
       
        
        for(j=0; j<B_num; j++){
            #ifdef DAlloc
            INT* H = calloc((proteinEntriesA[i].length+32)*(batches[j].lengths[31]+32)*32, sizeof(INT));
            #endif

            #ifdef SERIAL
                SWAVX_SeqToSeq_serial(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H, P, proteinEntriesA[i].length, proteinEntriesB[j].length, NumOfTest, (maxVal+i+A_num*j));
            #elif GPT
                SWAVX_GPT_Func(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H, P, proteinEntriesA[i].length, proteinEntriesB[j].length, (maxVal+i+A_num*j));
            #elif Query
            if(proteinEntriesA[i].length < batches[j].lengths[31])
                SWAVX_256_SeqToSeq_QueryLTBatch(batches[j], proteinEntriesA[i].protein, H, P, proteinEntriesA[i].length, NumOfTest, (maxVal+i+A_num*j), query_prof, (proteinEntriesA[i].length+2)*(batches[j].lengths[31]+2));
            else
                SWAVX_256_SeqToSeq_QueryBTBatch(proteinEntriesA[i].protein, batches[j], H, P, proteinEntriesA[i].length, NumOfTest, (maxVal+i+A_num*j), query_prof, (proteinEntriesA[i].length+2)*(batches[j].lengths[31]+2));
            #elif PARASAIL
                parasail_result_t *result = NULL;
                #ifndef QP
                #ifdef L8
                #ifdef P_scan
                    result = parasail_sw_scan_avx2_256_8(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, matrix);
                #elif P_diag
                    result = parasail_sw_diag_avx2_256_8(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, matrix);
                #else
                    result = parasail_sw_striped_avx2_256_8(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #endif
                #elif L16
                #ifdef P_scan
                    result = parasail_sw_scan_avx2_256_16(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #elif P_diag
                    result = parasail_sw_diag_avx2_256_16(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #else
                    result = parasail_sw_striped_avx2_256_16(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, &parasail_blosum62);
                #endif
                #else
                #ifdef P_scan
                    result = parasail_sw_scan_avx2_256_32(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, matrix);
                #elif P_diag
                    result = parasail_sw_diag_avx2_256_32(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, matrix);
                #else
                    result = parasail_sw_striped_avx2_256_32(proteinEntriesA[i].protein, proteinEntriesA[i].length, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1, matrix);
                #endif
                #endif
                #else
                    result = parasail_sw_scan_profile_avx2_256_8(profile, proteinEntriesB[j].protein, proteinEntriesB[j].length, 11, 1);
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
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j), query_prof,0);
                        else
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j), query_prof,1);
                        #else
                        if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j), query_prof,0);
                        else
                            SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j), query_prof,1);
                        #endif
                    #else
                    if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                        SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, (maxVal+i+A_num*j), query_prof,0);
                    else
                        SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, (maxVal+i+A_num*j), query_prof,1);
                    #endif
                #endif
            #endif
            #ifdef SAVEHP
            start += (proteinEntriesA[i].length+1) * (proteinEntriesB[j].length+1);
            #endif
            #ifdef DAlloc
            free(H);
            #endif
        }
        
        #ifdef PARASAIL
        #ifdef QP
        parasail_profile_free(profile);
        #endif
        #endif
    }
    
}