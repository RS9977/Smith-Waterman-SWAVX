#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>


#ifdef B512
#include "SWAVX_512_Func_SeqToSeq_SubMat.h"
#else
#include "SWAVX_256_Func_SeqToSeq_SubMat.h"
#endif

#define NumOfTest 1e0
//#define MaxHSize 3e7
#define bestMatchNum 2

typedef struct {
    ProteinEntry *proteinA;
    ProteinEntry *proteinB;
    INT* H;
    INT* P;
    int  A_num;
    int  B_num;
    INT* maxVal;
} WorkerIns;

void* chunck_computations(void* in);