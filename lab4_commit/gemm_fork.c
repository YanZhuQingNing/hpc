#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/wait.h>  // 包含wait头文件

#define N 1024  // 矩阵的大小

// 随机初始化矩阵
void initialize_matrix(double **matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            matrix[i][j] = (double)rand() / RAND_MAX * 10.0;  // 随机生成0到10之间的浮点数
        }
    }
}

// 简单的矩阵乘法函数
void naive_gemm(double **A, double **B, double **C, int start, int end, int size) {
    for (int i = start; i < end; i++) {
        for (int j = 0; j < size; j++) {
            C[i - start][j] = 0.0;  // 初始化C矩阵的元素
            for (int k = 0; k < size; k++) {
                C[i - start][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// 合并两个矩阵
void merge_matrices(double **C1, double **C2, double **C, int size) {
    for (int i = 0; i < size / 2; i++) {
        for (int j = 0; j < size; j++) {
            C[i][j] = C1[i][j];
        }
    }
    for (int i = size / 2; i < size; i++) {
        for (int j = 0; j < size; j++) {
            C[i][j] = C2[i - size / 2][j];
        }
    }
}

// 获取程序运行时间
void get_time(struct timeval *start, struct timeval *end, struct rusage *usage) {
    gettimeofday(end, NULL);
    getrusage(RUSAGE_SELF, usage);

    double real_time = (end->tv_sec - start->tv_sec) + (end->tv_usec - start->tv_usec) / 1000000.0;
    double user_time = usage->ru_utime.tv_sec + usage->ru_utime.tv_usec / 1000000.0;
    double sys_time = usage->ru_stime.tv_sec + usage->ru_stime.tv_usec / 1000000.0;

    printf("Real time: %f seconds\n", real_time);
    printf("User time: %f seconds\n", user_time);
    printf("Sys time: %f seconds\n", sys_time);
    printf("(sys + user) / real: %f\n", (user_time + sys_time) / real_time);
}

int main() {
    srand(time(NULL));  // 初始化随机数种子

    double **A = malloc(N * sizeof(double *));
    double **B = malloc(N * sizeof(double *));
    double **C = malloc(N * sizeof(double *));
    double **C1 = malloc((N / 2) * sizeof(double *));
    double **C2 = malloc((N / 2) * sizeof(double *));

    for (int i = 0; i < N; i++) {
        A[i] = malloc(N * sizeof(double));
        B[i] = malloc(N * sizeof(double));
        C[i] = malloc(N * sizeof(double));
        if (i < N / 2) {
            C1[i] = malloc(N * sizeof(double));
            C2[i] = malloc(N * sizeof(double));
        }
    }

    // 初始化矩阵
    initialize_matrix(A, N);
    initialize_matrix(B, N);
    initialize_matrix(C, N);

    struct timeval start, end;
    struct rusage usage;

    gettimeofday(&start, NULL);

    // 调用fork创建子进程
    pid_t pid = fork();

    if (pid < 0) {
        // fork失败
        fprintf(stderr, "Fork failed\n");
        return 1;
    } else if (pid == 0) {
        // 子进程
        naive_gemm(A, B, C1, 0, N / 2, N);
        exit(0);  // 子进程完成任务后退出
    } else {
        // 父进程
        naive_gemm(A, B, C2, N / 2, N, N);

        // 等待子进程完成
        wait(NULL);

        // 合并结果
        merge_matrices(C1, C2, C, N);

        gettimeofday(&end, NULL);
        getrusage(RUSAGE_SELF, &usage);

        // 打印时间
        get_time(&start, &end, &usage);
    }

    // 释放内存
    for (int i = 0; i < N; i++) {
        free(A[i]);
        free(B[i]);
        free(C[i]);
        if (i < N / 2) {
            free(C1[i]);
            free(C2[i]);
        }
    }
    free(A);
    free(B);
    free(C);
    free(C1);
    free(C2);

    return 0;
}
