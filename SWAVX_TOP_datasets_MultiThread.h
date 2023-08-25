#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>

#ifdef B512
#include "SWAVX_512_Func_SeqToSeq_SubMat.h"
#else
#include "SWAVX_256_Func_SeqToSeq_SubMat.h"
#endif

#define NumOfTest 1e0

typedef struct {
    ProteinEntry *proteinA;
    ProteinEntry *proteinB;
    int* H;
    int* P;
    int  A_num;
    int  B_num;
} WorkerIns;

void* chunck_computations(void* in);