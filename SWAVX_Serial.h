#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/stat.h>

#include "SWAVX_utils.h"
#include "SWAVX_SubMat.h"

void SWAVX_SeqToSeq_serial(int8_t *a, int8_t *b, INT *H, INT* P, int m, int n, int NumOfTest, INT* maxVal);
void similarityScore(int8_t* a, int8_t* b, int m, int n, int i, int j, INT* H, INT* P, int* maxVal);
int matchMissmatchScore(int i, int j, int8_t*a, int8_t*b);