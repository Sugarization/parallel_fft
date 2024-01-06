#include "fft_parallel.h"
#ifndef NDEBUG
#include <cstdio>
#endif

void directParallelFFT(NComplex *F, NComplex *x, indexT N, int inv)
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
    for (int k = 0; k < K; ++ k) {
        indexT m = 1ull << k;
        #pragma omp parallel for 
        for (indexT t = 0; t < N; ++ t) {
            indexT b = m << 1;
            indexT i = t & (b - 1);     // t mod 2m
            if (i < m) {  

                T phase = (T)inv * i * PI / m;
                T &R = x[t + m].re;
                T &I = x[t + m].im;
                NComplex v {R * cos(phase) - I * sin(phase), R * sin(phase) + I * cos(phase)};
                F[t] =     x[t] + v;
                F[t + m] = x[t] - v;
            }
        }
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