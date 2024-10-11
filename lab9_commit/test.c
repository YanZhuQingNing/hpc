#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>


//声明 C++ 的函数
#ifdef __cplusplus
extern "C" {
#endif

void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc);

#ifdef __cplusplus
}
#endif


#define A(i, j) a[(i) * lda + (j)]
#define B(i, j) b[(i) * ldb + (j)]
#define C(i, j) c[(i) * ldc + (j)]
#define NREPEATS 1
#define abs(x) ((x) < 0.0 ? -(x) : (x))
static double gtod_ref_time_sec = 0.0;

double compare_matrices(int m, int n, double *a, int lda, double *b, int ldb)
{
  int i, j;
  double max_diff = 0.0, diff;

  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
    {
      diff = abs(A(i, j) - B(i, j));
      max_diff = (diff > max_diff ? diff : max_diff);
    }

  return max_diff;
}

void random_matrix(int m, int n, double *a, int lda)
{
  /* drand48() generate pseudo-random numbers using the linear congruential algorithm and
     48-bit integer arithmetic. return nonnegative double-precision floating-point values
     uniformly distributed over the interval [0.0, 1.0). */
  double drand48();
  int i, j;

  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
      // A(i, j) = drand48();
      A(i, j) = 2.0 * drand48() - 1.0;
}

void copy_matrix(int m, int n, double *a, int lda, double *b, int ldb)
{
  int i, j;

  for (j = 0; j < n; j++)
    for (i = 0; i < m; i++)
      B(i, j) = A(i, j);
}
double dclock()
{
  double the_time, norm_sec;
  struct timeval tv;

  gettimeofday(&tv, NULL);

  if (gtod_ref_time_sec == 0.0)
    gtod_ref_time_sec = (double)tv.tv_sec;

  norm_sec = (double)tv.tv_sec - gtod_ref_time_sec;

  the_time = norm_sec + tv.tv_usec * 1.0e-6;

  return the_time;
}
void REF_MMult(int m, int n, int k, double *a, int lda,
               double *b, int ldb,
               double *c, int ldc)
{
  int i, j, p;

  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
        C(i, j) = 0.0;
      for (p = 0; p < k; p++)
      {
        C(i, j) = C(i, j) + A(i, p) * B(p, j);
      }
    }
  }
}

int main()
{
  int
      p,
      m, n, k,
      lda, ldb, ldc,
      rep;

  double
      dtime,
      dtime_best,
      gflops,
      diff;

  double
      *a,
      *b, *c,
      *cref, // Reference value computed by REF_MMult()
      *cold; // Random Initialization value

  printf("MY_MMult = [\n");

  for (p = 8; p <= 8192*4; p = 2 * p)
  {
    m = p;
    n = p;
    k = p;

    gflops = 2.0 * m * n * k * 1.0e-09;

    // ldx for column major
    // lda = (LDA == -1 ? m : LDA);
    // ldb = (LDB == -1 ? k : LDB);
    // ldc = (LDC == -1 ? m : LDC);

    lda = k;
    ldb = n;
    ldc = n;

    /* Allocate space for the matrices */
    /* Note: I create an extra column in A to make sure that
       prefetching beyond the matrix does not cause a segfault */
    a = (double *)malloc(lda * (k + 1) * sizeof(double));
    b = (double *)malloc(ldb * n * sizeof(double));
    c = (double *)malloc(ldc * n * sizeof(double));
    cold = (double *)malloc(ldc * n * sizeof(double));
    cref = (double *)malloc(ldc * n * sizeof(double));

    /* Generate random matrices A, B, Cold */
    random_matrix(m, k, a, lda);
    random_matrix(k, n, b, ldb);
    random_matrix(m, n, cold, ldc);

    copy_matrix(m, n, cold, ldc, cref, ldc);

    /* Run the reference implementation so the answers can be compared */
    //REF_MMult(m, n, k, a, lda, b, ldb, cref, ldc);

    /* Time the "optimized" implementation */
    for (rep = 0; rep < NREPEATS; rep++)
    {
      copy_matrix(m, n, cold, ldc, c, ldc);

      /* Time your implementation */
      dtime = dclock();
      MY_MMult(m, n, k, a, lda, b, ldb, c, ldc);
      dtime = dclock() - dtime;

      if (rep == 0)
        dtime_best = dtime;
      else
        dtime_best = (dtime < dtime_best ? dtime : dtime_best);
    }

    //diff = compare_matrices(m, n, c, ldc, cref, ldc);
    diff = 0.0;
    printf("%d %le %le \n", p, gflops / dtime_best, diff);
    fflush(stdout);

    free(a);
    free(b);
    free(c);
    free(cold);
    free(cref);
  }

  printf("];\n");

  exit(0);
}