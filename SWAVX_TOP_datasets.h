#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#ifdef B512
#include "SWAVX_512_Func_SeqToSeq_SubMat.h"
#else
#include "SWAVX_256_Func_SeqToSeq_SubMat.h"
#endif

#define NumOfTest 1e0
#define bestMatchNum 2