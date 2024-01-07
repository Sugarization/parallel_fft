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
        if (0) {
            #pragma omp parallel for 
            for (indexT t = 0; t < N; t += 4) {
                if (t & m) continue;        
                  
                indexT i = t & (m - 1);
                
                T c0 = cos((inv * (int64_t)i) * PI / m);
                T c1 = cos((inv * (int64_t)(i + 1)) * PI / m);
                T c2 = cos((inv * (int64_t)(i + 2)) * PI / m);
                T c3 = cos((inv * (int64_t)(i + 3)) * PI / m);
                T s0 = sin((inv * (int64_t)i) * PI / m);
                T s1 = sin((inv * (int64_t)(i + 1)) * PI / m);
                T s2 = sin((inv * (int64_t)(i + 2)) * PI / m);
                T s3 = sin((inv * (int64_t)(i + 3)) * PI / m);
                
                NComplex v0 (x[t + m + 0].re * c0 - x[t + m + 0].im * s0, x[t + m + 0].re * s0 + x[t + m + 0].im * c0);
                NComplex v1 (x[t + m + 1].re * c1 - x[t + m + 1].im * s1, x[t + m + 1].re * s1 + x[t + m + 1].im * c1);
                NComplex v2 (x[t + m + 2].re * c2 - x[t + m + 2].im * s2, x[t + m + 2].re * s2 + x[t + m + 2].im * c2);
                NComplex v3 (x[t + m + 3].re * c3 - x[t + m + 3].im * s3, x[t + m + 3].re * s3 + x[t + m + 3].im * c3);
                F[t + 0] = x[t + 0] + v0;
                F[t + 1] = x[t + 1] + v1;
                F[t + 2] = x[t + 2] + v2;
                F[t + 3] = x[t + 3] + v3;
                F[t + m + 0] = x[t + 0] - v0;
                F[t + m + 1] = x[t + 1] - v1;
                F[t + m + 2] = x[t + 2] - v2;
                F[t + m + 3] = x[t + 3] - v3;
            }
        } else {
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