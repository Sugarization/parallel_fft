#include "fft_parallel.h"
#ifndef NDEBUG
#include <cstdio>
#endif

void parallelBF1(NComplex *F, NComplex *x, indexT N)
{
    #pragma omp parallel for
    for (indexT t = 0; t < (N >> 1); ++ t) {
        indexT t0 = t << 1, t1 = (t << 1) | 1;
        F[t0] = x[t0] + x[t1];
        F[t1] = x[t0] - x[t1];
    }
}

void parallelBF2(NComplex *F, NComplex *x, indexT N, int inv)
{
    if (inv == -1) {
        #pragma omp parallel for
        for (indexT t = 0; t < (N >> 2); ++ t) {
            indexT t0 = (t << 2) | 0;
            indexT t1 = (t << 2) | 1;
            indexT t2 = (t << 2) | 2;
            indexT t3 = (t << 2) | 3;
            F[t0] = x[t0] + x[t2];
            F[t1] = x[t1] - x[t3].j();
            F[t2] = x[t0] - x[t2];
            F[t3] = x[t1] + x[t3].j();
        }
    } else {
        #pragma omp parallel for
        for (indexT t = 0; t < (N >> 2); ++ t) {
            indexT t0 = (t << 2) | 0;
            indexT t1 = (t << 2) | 1;
            indexT t2 = (t << 2) | 2;
            indexT t3 = (t << 2) | 3;
            F[t0] = x[t0] + x[t2];
            F[t1] = x[t1] + x[t3].j();
            F[t2] = x[t0] - x[t2];
            F[t3] = x[t1] - x[t3].j();
        }
    }
}

void parallelButterfly(NComplex *F, NComplex *x, indexT N, indexT k, int inv)
{
    if (k == 0) {
        parallelBF1(F, x, N);
        return;
    } else if (k == 1) {
        parallelBF2(F, x, N, inv);
        return;
    }
    indexT m = 1ull << k;
    #pragma omp parallel for 
    for (indexT r = 0; r < N / 2; r ++) {
        indexT t = ((r >> k) << (k + 1)) | r & (m - 1);
        indexT i = t & ((m << 1) - 1);     // t mod 2m
        if (i < m) { 
            NComplex v = x[t + m] * expJ((T)inv * i * PI / m);
            F[t] =     x[t] + v;
            F[t + m] = x[t] - v;
        }
    }
}

void parallelBitReversePermute(NComplex *x, indexT N)
{
    int K = MSB(N);
    #pragma omp parallel for 
    for (indexT i = 0; i < N; i ++) {
        indexT j = bitReverse(i, K);
        if (i < j) {
            NComplex tempor = x[i];
            x[i] = x[j];
            x[j] = tempor;
        }
    }
}

void directParallelFFT(NComplex *F, NComplex *x, indexT N, int inv)
{
    parallelBitReversePermute(x, N);

    int K = MSB(N);

    for (int k = 0; k < K; ++ k) {
        
        parallelButterfly(F, x, N, k, inv);
        
        NComplex *tmp = x;
        x = F;
        F = tmp;
    }
    if (!(K & 1)) {
        #pragma omp parallel for
        for (int i = 0; i < N; ++ i) {
            F[i] = x[i];
        }
    }
}

indexT bit2Reverse(indexT x_in, int K)
{
    return 0;
}

void embeddingParallelFFT(NComplex *F, NComplex *x, indexT N, int inv)
{
    int K = MSB(N);
    int K1 = K / 2;
    int K2 = K - K1;
    indexT N1 = (1lu << K1), N2 = (1lu << K2);

    #pragma omp parallel for 
    for (indexT j = 0; j < N2; ++ j) {
        iterativeColumnFFT(x, j, N1, N2, inv);
        for (indexT i = 0; i < N1; ++ i) {
            x[P2RC(i, j, K2)] = x[P2RC(i, j, K2)] * expJ(inv * 2.0 * PI * (i * j) / (T)N);  // twiddle factor
        }
    }

    #pragma omp parallel for 
    for (indexT i = 0; i < N1; ++ i) {
        iterativeRowFFT(x, i, N1, N2, inv);
    }
    parallelTranspose(F, x, N1, N2);
}

void embeddingParallelColumnCopyFFT(NComplex *F, NComplex *x, indexT N, int inv)
{
    int K = MSB(N);
    int K1 = K / 2;
    int K2 = K - K1;
    indexT N1 = (1lu << K1), N2 = (1lu << K2);

    #pragma omp parallel for 
    for (indexT j = 0; j < N2; ++ j) {
        for (indexT i = 0; i < N1; ++ i) {
            F[P2RC(j, i, K1)] = x[P2RC(i, j, K2)];
        }
        iterativeRowFFT(F, j, N2, N1, inv);
        for (indexT i = 0; i < N1; ++ i) {
            F[P2RC(j, i, K1)] = F[P2RC(j, i, K1)] * expJ(inv * 2.0 * PI * (i * j) / (T)N);  // twiddle factor
        }
    }
    parallelTranspose(x, F, N2, N1);

    #pragma omp parallel for 
    for (indexT i = 0; i < N1; ++ i) {
        iterativeRowFFT(x, i, N1, N2, inv);
    }
    parallelTranspose(F, x, N1, N2);
}