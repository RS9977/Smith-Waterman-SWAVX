#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3

#define RESET   "\033[0m"
#define BOLDRED "\033[1m\033[31m"      /* Bold Red */

#define BUFFER_SIZE 128
#define MAX_THREAD 40

typedef struct {
    int8_t *protein;
    int length;
} ProteinEntry;

#ifdef L8
typedef int8_t INT;
#elif L16
typedef short int INT;
#else
typedef int INT;
#endif



void printMatrix(INT* matrix, int8_t *a, int8_t *b, int m, int n);
void printPredecessorMatrix(INT* matrix, int8_t *a, int8_t *b, int m, int n);
void generate(int8_t *a, int8_t *b, int m, int n);
void saveInFile (INT *H, int8_t *a, int8_t *b, int m, int n);
int readProteinDataset(const char *filename, ProteinEntry **proteinEntries, int *numEntries);
int getNumCPUThreads();
void load_balance(int* chunck_start, int* chunck_num, int* chunck_size, int Hsize, int numEntries, ProteinEntry *proteinEntries, int NumOfThreads);
double wakeup_delay();
double interval(struct timespec start, struct timespec end);
