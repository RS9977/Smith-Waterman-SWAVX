#include "SWAVX_TOP_datasets.h"

int main(int argc, char* argv[]) {
    
    const char *filenameA = "test2.fasta"; 
    const char *filenameB = "test3.fasta"; 

    if(argc>1){
        filenameA = argv[1];
        filenameB = argv[2];
    }

    
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

    for(i=0; i<numEntriesA; i++)
        HsizeA += proteinEntriesA[i].length;
    for(i=0; i<numEntriesB; i++)
        HsizeB += proteinEntriesB[i].length;

    //Allocates similarity matrix H
    INT *H;
    H = calloc((HsizeA+numEntriesA) * (HsizeB+numEntriesB), sizeof(INT));

    //Allocates predecessor matrix P
    INT *P;
    #ifdef BT
    P = calloc((HsizeA+numEntriesA) * (HsizeB+numEntriesB), sizeof(INT));
    #endif

    struct timespec time_start, time_stop;
    clock_gettime(CLOCK_REALTIME, &time_start);

    long long int start = 0;
    for(i=0; i<numEntriesA; i++){
        for(j=0; j<numEntriesB; j++){
            #ifdef B512
                #ifdef BT
                if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                    SWAVX_512_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, gapscore);
                else
                    SWAVX_512_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, gapscore);
                #else
                if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                    SWAVX_512_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, gapscore);
                else
                    SWAVX_512_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, gapscore);
                #endif
            #else
                #ifdef BT
                if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                    SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P+start, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, gapscore);
                else
                    SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P+start, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, gapscore);
                #else
                if(proteinEntriesA[i].length > proteinEntriesB[j].length)
                    SWAVX_256_SeqToSeq_SubMat(proteinEntriesA[i].protein, proteinEntriesB[j].protein, H+start, P, proteinEntriesA[i].length+1, proteinEntriesB[j].length+1, NumOfTest, gapscore);
                else
                    SWAVX_256_SeqToSeq_SubMat(proteinEntriesB[j].protein, proteinEntriesA[i].protein, H+start, P, proteinEntriesB[j].length+1, proteinEntriesA[i].length+1, NumOfTest, gapscore);
                #endif
            #endif
            start += (proteinEntriesA[i].length+1) * (proteinEntriesB[j].length+1);
        }
    }

    //Gets final time
    clock_gettime(CLOCK_REALTIME, &time_stop);
    double MeanTime = interval(time_start, time_stop)/NumOfTest;
    printf("Elapsed time: %f\n", MeanTime);
    printf("GCUPS: %f\n", HsizeA*HsizeB/(1e9*MeanTime));


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
