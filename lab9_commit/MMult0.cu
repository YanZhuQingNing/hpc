#include <stdio.h>

// CUDA 核函数进行矩阵乘法
__global__ void matrixMultiply(double *A, double *B, double *C, int m, int n, int k) {
    int row = blockIdx.y * blockDim.y + threadIdx.y; // 当前处理的行
    int col = blockIdx.x * blockDim.x + threadIdx.x; // 当前处理的列

    if (row < m && col < n) {
        double sum = 0.0;
        for (int p = 0; p < k; p++) {
            sum += A[row * k + p] * B[p * n + col];
        }
        C[row * n + col] = sum;
    }
}

// 主机函数，调用 CUDA 核函数
extern "C" void MY_MMult(int m, int n, int k, double *a, int lda,
                         double *b, int ldb,
                         double *c, int ldc) {
    double *d_A, *d_B, *d_C;

    size_t sizeA = m * k * sizeof(double);
    size_t sizeB = k * n * sizeof(double);
    size_t sizeC = m * n * sizeof(double);

    // 在设备上分配内存
    cudaMalloc(&d_A, sizeA);
    cudaMalloc(&d_B, sizeB);
    cudaMalloc(&d_C, sizeC);

    // 将 A 和 B 从主机复制到设备
    cudaMemcpy(d_A, a, sizeA, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, b, sizeB, cudaMemcpyHostToDevice);

    // 定义线程块和网格的维度
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((n + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (m + threadsPerBlock.y - 1) / threadsPerBlock.y);

    // 调用 CUDA 核函数
    matrixMultiply<<<numBlocks, threadsPerBlock>>>(d_A, d_B, d_C, m, n, k);

    // 将结果从设备复制回主机
    cudaMemcpy(c, d_C, sizeC, cudaMemcpyDeviceToHost);

    // 释放设备内存
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}
