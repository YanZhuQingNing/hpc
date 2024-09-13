#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc) {
    #pragma omp parallel for
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int p = 0; p < k; p++) {
                c[i*ldc+j] += a[i * lda + p] * b[p * ldb + j];
            }
        }
    }
}
