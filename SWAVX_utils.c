#include "SWAVX_utils.h"

/*--------------------------------------------------------------------
 * Function:    printMatrix
 * Purpose:     Print Matrix
 */
void printMatrix(int* matrix, int8_t *a, int8_t *b, int m, int n) {
    long long int i, j, ind;
    printf(" \t \t");
    for(i=0; i<m-1; i++)
        printf("%d\t",a[i]);
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                printf(" \t");
            else if(j==-1 && i>0)
                printf("%d\t",b[i-1]); 
            else
                printf("%d\t", matrix[ind]);
        }
        printf("\n");
    }
}  /* End of printMatrix */

/*--------------------------------------------------------------------
 * Function:    printPredecessorMatrix
 * Purpose:     Print predecessor matrix
 */
void printPredecessorMatrix(int* matrix, int8_t *a, int8_t *b, int m, int n) {
    long long int i, j, ind;
    printf("    ");
    for(i=0; i<m-1; i++)
        printf("%d ",a[i]);
    printf("\n");
    for (i = 0; i < n; i++) { //Lines
        for (j = -1; j < m; j++) {
            if(i+j<n)
                ind = (i+j)*(i+j+1)/2 + i;
            else if(i+j<m)
                ind = (n+1)*(n)/2 + (i+j-n)*n + i;
            else
                ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
            if(i+j<0)
                printf("  ");
            else if(j==-1 && i>0)
                printf("%d ",b[i-1]); 
            else{
                if(matrix[ind] < 0) {
                    printf(BOLDRED);
                    if (matrix[ind] == -UP)
                        printf("↑ ");
                    else if (matrix[ind] == -LEFT)
                        printf("← ");
                    else if (matrix[ind] == -DIAGONAL)
                        printf("↖ ");
                    else
                        printf("- ");
                    printf(RESET);
                } else {
                    if (matrix[ind] == UP)
                        printf("↑ ");
                    else if (matrix[ind] == LEFT)
                        printf("← ");
                    else if (matrix[ind] == DIAGONAL)
                        printf("↖ ");
                    else
                        printf("- ");
            }
            }
        }
        printf("\n");
    }

}  /* End of printPredecessorMatrix */


/*--------------------------------------------------------------------
 * Function:    generate
 * Purpose:     Generate arrays a and b
 */
 void generate(int8_t *a, int8_t *b, int m, int n){
    //Generates the values of a
    long long int i;
    for(i=0;i<m;i++){
        a[i] = rand()%4;
        /*
        int aux=rand()%4;
        if(aux==0)
            a[i]='A';
        else if(aux==2)
            a[i]='C';
        else if(aux==3)
            a[i]='G';
        else
            a[i]='T';
        */
    }

    //Generates the values of b
    for(i=0;i<n;i++){
        b[i] = rand()%4;
        /*
        int aux=rand()%4;
        if(aux==0)
            b[i]='A';
        else if(aux==2)
            b[i]='C';
        else if(aux==3)
            b[i]='G';
        else
            b[i]='T';
        */
    }
} /* End of generate */


//Saving the Similarity matrix in a file
void saveInFile (int *H, int8_t *a, int8_t *b, int m, int n){
    FILE* file = fopen("256_ST_SB.txt", "w");

    // Check if the file was opened successfully
    if (file == NULL) {
        printf("Failed to open the file.\n");
    }
    else{
        int ind;
        fprintf(file, " \t \t");
        int i, j;
        for(i=0; i<m-1; i++)
            fprintf(file,"%d\t",a[i]);
        fprintf(file, "\n");
        for (i = 0; i < n; i++) { //Lines
            for (j = -1; j < m; j++) {
                if(i+j<n)
                    ind = (i+j)*(i+j+1)/2 + i;
                else if(i+j<m)
                    ind = (n+1)*(n)/2 + (i+j-n)*n + i;
                else
                    ind = (i*j) + ((m-j)*(i+(i-(m-j-1))))/2 + ((n-i)*(j+(j-(n-i-1))))/2 + (m-j-1);
                if(i+j<0)
                    fprintf(file, " \t");
                else if(j==-1 && i>0)
                    fprintf(file, "%d\t",b[i-1]); 
                else
                    fprintf(file, "%d\t", H[ind]);
            }
            fprintf(file, "\n");
        }
    }
}


//Read Dataset

int readProteinDataset(const char *filename, ProteinEntry **proteinEntries, int *numEntries) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return 0;
    }

    int maxEntries = 1000; // Initial capacity, can be adjusted
    *proteinEntries = (ProteinEntry *)malloc(maxEntries * sizeof(ProteinEntry));
    if (*proteinEntries == NULL) {
        perror("Memory allocation failed");
        fclose(file);
        return 0;
    }

    int entryCount = 0;
    char line[1024]; // Adjust this based on your dataset
    int8_t *sequence = NULL;
    int seqLength = 0;
    int sequenceStarted = 0;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '>') {
            sequenceStarted = 1;
            if (sequence != NULL) {
                (*proteinEntries)[entryCount].protein = sequence;
                (*proteinEntries)[entryCount].length = seqLength;
                entryCount++;
            }
            sequence = NULL;
            seqLength = 0;
        } else if (sequenceStarted) {
            line[strcspn(line, "\r\n")] = 0; // Remove trailing newline characters
            int lineLength = strlen(line);
            if (lineLength > 0) {
                if (sequence == NULL) {
                    sequence = (int8_t *)malloc(lineLength * sizeof(int8_t)+8);
                } else {
                    sequence = (int8_t *)realloc(sequence, (seqLength + lineLength) * sizeof(int8_t)+8);
                }
                if (sequence == NULL) {
                    perror("Memory allocation failed");
                    fclose(file);
                    return 0;
                }
                for (int i = 0; i < lineLength; i++) {
                    sequence[seqLength + i] = line[i] - 'A'; // Convert char to numeric value
                }
                seqLength += lineLength;
            }
        }
    }

    if (sequence != NULL) {
        (*proteinEntries)[entryCount].protein = sequence;
        (*proteinEntries)[entryCount].length = seqLength;
        entryCount++;
    }

    fclose(file);
    *numEntries = entryCount;
    return 1;
}

//Get CPU threads
int getNumCPUThreads() {
    char buffer[BUFFER_SIZE];
    char* command = "lscpu | grep 'CPU(s):'";

    FILE* fp = popen(command, "r");
    if (fp == NULL) {
        perror("popen");
        return 0;
    }

    int numThreads=-1;
    while (fgets(buffer, BUFFER_SIZE, fp) != NULL) {
        sscanf(buffer, "CPU(s): %d", &numThreads);
    }

    pclose(fp);
    
    return numThreads;
}


//Load balance of dataset for each threads
void load_balance(int* chunck_start, int* chunck_num, int* chunck_size, int Hsize, int numEntries, ProteinEntry *proteinEntries, int NumOfThreads){
    int cnt = 0;
    int sizeOfChunk = Hsize/NumOfThreads;
    int temp_size = 0;
    int temp_temp_size = 0;
    chunck_start[0] = 0;
    int temp_i=-1;
    int i;
    for(i=0; i<numEntries; i++){
        temp_size += proteinEntries[i].length;
        if(temp_size>sizeOfChunk*(cnt+1)){
            chunck_start[cnt+1] = i+1;
            chunck_size[cnt] = temp_size - temp_temp_size;
            chunck_num[cnt]  = i - temp_i;
            cnt ++;
            temp_temp_size = temp_size;
            temp_i = i;
        }
        if(cnt==NumOfThreads-1){
            chunck_size[cnt] = Hsize - temp_temp_size;
            chunck_num[cnt]  = numEntries-1 -temp_i;
            break;
        }
    }
}

//Wake Up
double interval(struct timespec start, struct timespec end)
{
  struct timespec temp;
  temp.tv_sec = end.tv_sec - start.tv_sec;
  temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (temp.tv_nsec < 0) {
    temp.tv_sec = temp.tv_sec - 1;
    temp.tv_nsec = temp.tv_nsec + 1000000000;
  }
  return (((double)temp.tv_sec) + ((double)temp.tv_nsec)*1.0e-9);
}

double wakeup_delay()
{
  double meas = 0; int j, i;
  struct timespec time_start, time_stop;
  double quasi_random = 0;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_start);
  j = 1e7;
  while (meas < 1.0) {
    for (i=1; i<j; i++) {
      /* This iterative calculation uses a chaotic map function, specifically
         the complex quadratic map (as in Julia and Mandelbrot sets), which is
         unpredictable enough to prevent compiler optimisation. */
      quasi_random = quasi_random*quasi_random - 1.923432;
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_stop);
    meas = interval(time_start, time_stop);
    j *= 2; /* Twice as much delay next time, until we've taken 1 second */
  }
  return quasi_random;
}