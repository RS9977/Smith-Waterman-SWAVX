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
        NumOfThreads = strtoll(argv[1], NULL, 10);
    }

    //Reading datasets
    ProteinEntry *proteinEntriesA;
    int numEntriesA;
    if (!readProteinDataset(filenameA, &proteinEntriesA, &numEntriesA)) {
        printf("Unable to read the file: %s\n",filenameA);
        return 0;        
    }
    

    ProteinEntry *proteinEntriesB;
    int numEntriesB;
    if (!readProteinDataset(filenameB, &proteinEntriesB, &numEntriesB)) {
        printf("Unable to read the file: %s\n",filenameB);
        return 0;        
    }
    
    long long int HsizeA = 0;
    long long int HsizeB = 0;
    int i,j;

    for(i=0; i<numEntriesA; i++)
        HsizeA += proteinEntriesA[i].length;
    for(i=0; i<numEntriesB; i++)
        HsizeB += proteinEntriesB[i].length;
    
    //Load balancing
    int B_chunck_start[MAX_THREAD];
    int B_chunck_num  [MAX_THREAD];
    int B_chunck_size [MAX_THREAD];
    load_balance(B_chunck_start, B_chunck_num, B_chunck_size, HsizeB, numEntriesB, proteinEntriesB, NumOfThreads);

    //Allocates similarity matrix H
    int *H;
    H = calloc((HsizeA+numEntriesA) * (HsizeB+numEntriesB), sizeof(int));

    //Allocates predecessor matrix P
    int *P;
    P = calloc((HsizeA+numEntriesA) * (HsizeB+numEntriesB), sizeof(int));

    double initialTime = omp_get_wtime();

        pthread_t threads[MAX_THREAD];
        WorkerIns thread_data_array[MAX_THREAD];
        int t;
        int rc;          
        long long int start = 0;
        for(t=0; t<NumOfThreads; t++){ 
            
            thread_data_array[t].proteinA = proteinEntriesA;  
            thread_data_array[t].proteinB = proteinEntriesB + B_chunck_start[t];  
            thread_data_array[t].H        = H+start;
            thread_data_array[t].P        = P+start;
            thread_data_array[t].A_num    = numEntriesA;
            thread_data_array[t].B_num    = B_chunck_num[t];
                
            rc = pthread_create(&threads[t], NULL, chunck_computations, (void*) &thread_data_array[t]);
            if (rc) {
                printf("ERROR; return code from pthread_create() is %d\n", rc);
                exit(-1);
            }
            start += (HsizeA+numEntriesA) * (B_chunck_size[t]+B_chunck_num[t]);
        }
        
        for (t = 0; t<NumOfThreads; t++) {
            if (pthread_join(threads[t], NULL)){ 
                printf("ERROR; code on return from join is %d\n", rc);
                exit(-1);
            }
        }


    //Gets final time
    double finalTime = omp_get_wtime();

    double MeanTime = (finalTime - initialTime)/NumOfTest;
    printf("Elapsed time: %f\n", MeanTime);
    printf("GCUPS: %f\n", HsizeA*HsizeB/(1e9*MeanTime));
    printf("WakeUpVal: %f\n",wk);


    //Frees similarity matrixes
    free(H);
    free(P);


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
    int* H                        = inss -> H;
    int* P                        = inss -> P;
    int  A_num                    = inss -> A_num;
    int  B_num                    = inss -> B_num;
    long long int start           = 0;
    int  i,j;
    for(i=0; i<A_num; i++){
        for(j=0; j<B_num; j++){
            #ifdef B512
            if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                SWAVX_512_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, gapscore);
            else
                SWAVX_512_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, gapscore);
            #else
            if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, gapscore);
            else
                SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, gapscore);
            #endif
            start += (proteinEntriesA[i].length+1) * (proteinEntriesB[j].length+1);
        }
    }
}