#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define N 8  // 矩阵大小

void print_matrix(double *matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i * n + j]);
        }
        printf("\n");
    }
}
double compare_matrices(int m, int n, double *a, double *b)
{
  int i, j;
  double max_diff = 0.0, diff;

  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
    {
      diff = abs(a[i*m+j] - b[i*m+j]);
      max_diff = (diff > max_diff ? diff : max_diff);
    }

  return max_diff;
}
int main(int argc, char **argv) {
    int rank, size;
    double *A, *B, *C, *local_A, *local_C,*C_REF;
    int local_n, start_row, end_row;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 计算每个进程负责的行数
    local_n = N / size;
    start_row = rank * local_n;
    end_row = start_row + local_n;

    // 分配本地矩阵内存
    local_A = (double *)malloc(local_n * N * sizeof(double));
    local_C = (double *)malloc(local_n * N * sizeof(double));

    // 主进程初始化全局矩阵 A 和 B
    if (rank == 0) {
        A = (double *)malloc(N * N * sizeof(double));
        B = (double *)malloc(N * N * sizeof(double));
        C = (double *)malloc(N * N * sizeof(double));
        C_REF = (double *)malloc(N * N * sizeof(double));
        // 初始化矩阵 A 和 B
        // for (int i = 0; i < N * N; i++) {
        //     A[i] = i + 1;
        //     B[i] = i + 1;
        // }
        if(N==4){
            double A_data[]={
            -1.000000e+00, -2.707955e-01, 5.350056e-02, 8.634630e-01, 
            -9.980292e-01, -8.173388e-01, -9.113315e-02, 1.361192e-01, 
            -9.167380e-01, -8.154047e-01, -5.336431e-01, 1.121887e-01, 
            -6.467147e-01, -2.556555e-02, 6.625836e-01, -8.983362e-01, 
            };

            double B_data[]={
            5.341023e-01, 7.519617e-01, 6.208589e-01, -8.464509e-01, 
            -9.621704e-01, 6.311373e-02, -6.231595e-01, 6.305478e-01, 
            -4.952805e-01, 8.405219e-01, 7.726289e-01, 9.697820e-01, 
            -4.036057e-01, 3.086230e-02, 1.412280e-01, -7.632966e-01, 
            };

            double C_REF_data[]={
            -6.485472e-01, -6.974358e-01, -2.888289e-01, 6.850692e-02, 
            2.435674e-01, -8.744635e-01, -1.614912e-01, 1.371329e-01, 
            5.139494e-01, -1.185891e+00, -4.575017e-01, -3.413287e-01, 
            -2.864046e-01, 4.127302e-02, -5.262187e-04, 1.859551e+00, 
            };
            for (int i = 0; i < N * N; i++) {
                A[i] = A_data[i];
                B[i] = B_data[i];
                C_REF[i] = C_REF_data[i];
            }
        }
        else if(N==8){
            double A_data[]={
            1.900078e-02, 8.548636e-02, -3.314322e-01, 2.440604e-02, -6.461445e-01, -7.460251e-01, 5.706861e-01, -6.391164e-02, 
            -3.681220e-01, 3.359572e-01, -8.358069e-01, 1.730207e-01, -7.762620e-01, -9.918196e-01, -1.155304e-01, -6.060053e-01, 
            4.034051e-01, -6.577163e-01, -4.593707e-01, 4.775171e-01, 7.053617e-01, -9.989944e-01, -4.241922e-01, -1.056461e-01, 
            -2.655067e-01, 6.615583e-01, -8.514753e-01, 7.433776e-01, 6.223437e-01, -6.140962e-01, 7.067059e-01, 4.568796e-01, 
            7.526247e-01, 7.238108e-01, 3.524515e-01, 1.559346e-01, 9.097746e-01, 1.735603e-01, 5.194394e-02, 9.321860e-01, 
            -2.580107e-01, 5.384423e-01, 6.203171e-01, 9.294072e-01, 8.085664e-01, -2.500817e-01, 4.568334e-01, 5.136816e-01, 
            5.050480e-02, -8.948411e-01, -7.157754e-02, -3.329414e-01, -7.221695e-01, -9.464281e-01, -4.673714e-01, -6.788623e-02, 
            -1.443309e-01, 9.334593e-01, 1.383701e-01, -2.143441e-01, 3.723747e-01, -7.046941e-02, 5.669680e-01, 7.379329e-01, 
            };

            double B_data[]={
            -7.319215e-01, 8.165651e-01, -4.029224e-02, 8.852518e-01, 5.361704e-01, -2.435351e-01, -5.578819e-01, -3.161686e-01, 
            6.921052e-01, 5.488294e-01, -9.377434e-01, -2.239507e-01, 5.038825e-01, -1.625671e-01, 5.086909e-02, -8.637625e-02, 
            -9.638247e-01, -2.678706e-01, 2.168604e-01, 8.610245e-01, -8.182718e-01, -7.231037e-01, -1.094956e-02, 1.762435e-01, 
            8.817096e-01, 8.885760e-01, 5.845223e-01, 8.227684e-01, 7.079334e-01, -2.383066e-01, 7.650440e-01, 1.702999e-01, 
            3.401642e-01, 7.696433e-01, -6.350064e-01, -3.245633e-01, -6.507111e-01, -4.004399e-01, 5.476602e-01, -8.601459e-01, 
            5.999431e-01, -1.096578e-01, 5.024259e-02, 6.562465e-02, -9.604593e-01, -9.579758e-01, -8.819950e-01, 6.378629e-01, 
            -7.135664e-01, -6.739431e-01, 7.975177e-01, 4.041545e-01, -1.618930e-01, -5.223948e-01, -7.603972e-01, -1.936786e-01, 
            2.209206e-01, -1.503654e-01, -3.824202e-01, 9.420213e-01, -4.072602e-01, -9.257813e-01, 9.560251e-02, 1.709661e-01, 
            };

            double C_REF_data[]={
            -7.024897e-01, -6.175928e-01, 7.138586e-01, 6.358145e-02, 1.412361e+00, 9.497799e-01, -1.198863e-01, -1.091874e-01, 
            5.495471e-01, -5.828324e-02, 2.023827e-01, -1.409117e+00, 2.501541e+00, 2.480544e+00, 8.435441e-01, -7.664723e-02, 
            3.326345e-02, 1.469986e+00, -1.598580e-02, -6.368559e-02, 1.211029e+00, 1.221015e+00, 1.691703e+00, -1.250213e+00, 
            1.568243e+00, 1.036256e+00, -3.969600e-01, -3.099133e-02, 1.298360e+00, -5.740635e-02, 1.148578e+00, -9.824463e-01, 
            3.303513e-01, 1.561958e+00, -1.425547e+00, 1.551175e+00, -5.565087e-01, -2.013690e+00, 1.271708e-01, -7.343229e-01, 
            6.956050e-01, 1.009126e+00, -1.748624e-01, 1.339494e+00, -2.857656e-01, -1.493150e+00, 1.240700e+00, -5.529886e-01, 
            -1.375815e+00, -8.533851e-01, 6.912177e-01, -1.710138e-01, 8.812948e-01, 1.767113e+00, 4.605126e-01, 8.840335e-02, 
            2.721833e-01, -3.181246e-02, -1.034846e+00, 4.047740e-01, -4.389420e-01, -1.226528e+00, -1.319781e-01, -3.960062e-01, 
            };
            for (int i = 0; i < N * N; i++) {
                A[i] = A_data[i];
                B[i] = B_data[i];
                C_REF[i] = C_REF_data[i];
            }
        }
        
        
    }
    
        // 分发矩阵 A 的数据到各个进程
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            int start_row = p * local_n;
            int end_row = start_row + local_n;
            MPI_Send(&A[start_row * N], local_n * N, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
        }
        // 主进程复制自己的部分
        for (int i = 0; i < local_n; i++) {
            for (int j = 0; j < N; j++) {
                local_A[i * N + j] = A[i * N + j];
            }
        }
    } else {
        MPI_Recv(local_A, local_n * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

   if (rank == 0) {
        MPI_Bcast(B, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        B = (double *)malloc(N * N * sizeof(double));
        MPI_Bcast(B, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    // 每个进程计算自己的矩阵块
    for (int i = 0; i < local_n; i++) {
        for (int j = 0; j < N; j++) {
            local_C[i * N + j] = 0.0;
            for (int k = 0; k < N; k++) {
                local_C[i * N + j] += local_A[i * N + k] * B[k * N + j];
            }
        }
    }

    // 主进程收集计算结果
    if (rank == 0) {
        // 将主进程的结果复制到结果矩阵
        for (int i = 0; i < local_n; i++) {
            for (int j = 0; j < N; j++) {
                C[i * N + j] = local_C[i * N + j];
            }
        }

        // 接收其他进程的结果
        for (int p = 1; p < size; p++) {
            MPI_Recv(local_C, local_n * N, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < local_n; i++) {
                for (int j = 0; j < N; j++) {
                    C[(start_row + local_n * p + i) * N + j] = local_C[i * N + j];
                }
            }
        }
        double diff=compare_matrices(N,N,C, C_REF);
        printf("diff=%e\n",diff);
        //输出结果
        printf("Matrix A:\n");
        print_matrix(A, N);
        printf("Matrix B:\n");
        print_matrix(B, N);
        printf("Result Matrix C:\n");
        print_matrix(C, N);
        free(A);
        free(B);
        free(C);
    } else {
        // 其他进程发送计算结果给主进程
        MPI_Send(local_C, local_n * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    free(local_A);
    free(local_C);

    MPI_Finalize();
    return 0;
}
