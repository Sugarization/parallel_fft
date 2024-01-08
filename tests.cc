#include "tests.h"

void fft_test(NComplex *F, NComplex *x, indexT N, FFT_Type type, int nThreads, int inv)
{
    if (nThreads > 0) {
        omp_set_num_threads(nThreads);
    }
    
    srand(time(nullptr));
   
    
    // for (indexT i = 0; i < N; ++ i) {
        
    //     T error = (F1[i] - F[i]).abs() / x1[i].abs();
    //     if (error > 1e-6 && x1[i].abs() > 1e-6) {
    //         printf("Error check @ %d failed: %.8f, from %f, %f; %f, %f\n", (int) i, error, F1[i].re, F1[i].im, F[i].re, F[i].im);
    //     }
    // }
}

void transpose_test()
{
    NComplex x[30000], y[30000], z[30000];
    int N = 128;
    for (int i = 0; i < N; ++ i) for (int j = 0; j < N; ++ j) x[i*N + j] = NComplex(i*i, j);
    parallelTranspose(y, x, N, N);
    parallelNaiveTranspose(z, x, N, N);
    for (int i = 0; i < N; ++ i) {for (int j = 0; j < N; ++ j) if ((z[i*N+j] - y[i*N+j]).abs() > 1e-9) puts("err");}
}