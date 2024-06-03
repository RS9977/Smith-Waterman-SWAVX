#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>


#ifdef PARASAIL
#include "parasail.h"
#include "parasail/matrices/blosum62.h"
#endif


#ifdef SERIAL
#include "SWAVX_Serial.h"
#elif GPT
#include "SWAVX_GPT.h"
#elif GPTOPT
#include "SWAVX_GPT_Opt.h"
#elif B512
#include "SWAVX_512_Func_SeqToSeq_SubMat.h"
#elif Query
#include "SWAVX_256_Func_SeqToSeq_Query.h"
#else
#include "SWAVX_256_Func_SeqToSeq_SubMat.h"
#endif

#define NumOfTest 1e1
//#define MaxHSize 3e7
#define bestMatchNum 2
#define MaxLenQuery 2000
typedef struct {
    ProteinEntry *proteinA;
    #ifndef Query
    ProteinEntry *proteinB;
    #else
    ProteinBatch* batch;
    int MaxHSize;
    #endif
    #ifndef DAlloc
    INT* H;
    INT* P;
    #endif
    int  A_num;
    int  B_num;
    INT* maxVal;
} WorkerIns;

void* chunck_computations(void* in);