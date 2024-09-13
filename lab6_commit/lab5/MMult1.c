#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

typedef struct {
    int startRow;
    int endRow;
    int ROWS;
    int COLS;
    int colsA;
    double *matrixA;
    double *matrixB;
    double *matrixC;
} MatrixThreadArgs;

void *matrixMultiply(void *args) {
    MatrixThreadArgs *threadArgs = (MatrixThreadArgs *)args;
    int startRow = threadArgs->startRow;
    int endRow = threadArgs->endRow;
    int COLS = threadArgs->COLS;
    int colsA = threadArgs->colsA;
    double *matrixA = threadArgs->matrixA;
    double *matrixB = threadArgs->matrixB;
    double *matrixC = threadArgs->matrixC;

    for (int i = startRow; i < endRow; i++) {
        for (int j = 0; j < COLS; j++) {
            for (int k = 0; k < colsA; k++) {
                matrixC[i * COLS + j] += matrixA[i * colsA + k] * matrixB[k * COLS + j];
            }
        }
    }

    free(threadArgs);
    pthread_exit(NULL);
}

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc) {
    fprintf(stderr, __FILE__ "\n");
    pthread_t threads[8];
    MatrixThreadArgs *threadArgs[8];
    int rowsPerThread = m / 8;

    for (int thread_id = 0; thread_id < 8; thread_id++) {
        int startRow = thread_id * rowsPerThread;
        int endRow = (thread_id == 7) ? m : startRow + rowsPerThread; // 处理最后一个线程的边界情况
        threadArgs[thread_id] = (MatrixThreadArgs *)malloc(sizeof(MatrixThreadArgs));
        threadArgs[thread_id]->startRow = startRow;
        threadArgs[thread_id]->endRow = endRow;
        threadArgs[thread_id]->ROWS = m;
        threadArgs[thread_id]->COLS = n;
        threadArgs[thread_id]->colsA = k;
        threadArgs[thread_id]->matrixA = a;
        threadArgs[thread_id]->matrixB = b;
        threadArgs[thread_id]->matrixC = c;
        pthread_create(&threads[thread_id], NULL, matrixMultiply, threadArgs[thread_id]);
    }

    for (int i = 0; i < 8; i++) {
        if (pthread_join(threads[i], NULL) != 0) {
            perror("pthread_join");
            exit(EXIT_FAILURE);
        }
    }
}
