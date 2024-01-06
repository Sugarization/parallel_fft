#include "fft_base.h"
#include "fft_parallel.h"
#include <cstdlib>
#include <ctime>
#include <cstdio>

const indexT N = 1 << 25;
const T Fs = 64;

double difftime(timespec &b, timespec &a) {
    uint64_t diff = (b.tv_sec - a.tv_sec) * 1000000000ull + (b.tv_nsec - a.tv_nsec);
    return diff / 1e9;
}

int main()
{
    NComplex *x = new NComplex[N];
    NComplex *F = new NComplex[N];
    
    srand(time(nullptr));
    T f = 10;
    for (indexT i = 0; i < N; ++ i) {
        // x[i].re = sin(2 * PI * f / Fs * i);
        x[i].re = i;
        x[i].im = 0;
    }
    #ifndef NDEBUG
    NComplex *x1 = new NComplex[N];
    for (indexT i = 0; i < N; ++ i) {
        x1[i] = x[i];
        // printf("%.4f %.4f\n", x1[i].re, x1[i].im);
    }
    #endif
    // naiveDFT(F, x, N);
    timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    directParallelFFT(F, x, N);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    printf("Time elapsed: %.6f\n", difftime(t1, t0));

    #ifndef NDEBUG
    
    iterativeFFT(x1, N);
    for (indexT i = 0; i < N; ++ i) {
        
        T error = (x1[i] - F[i]).abs() / x1[i].abs();
        if (error > 1e-6 && x1[i].abs() > 1e-6) {
            printf("Error check @ %d failed: %.8f, from %f, %f; %f, %f\n", (int) i, error, x1[i].re, x1[i].im, F[i].re, F[i].im);
        }
    }
    delete[] x1;
    #endif

    delete[] x;
    delete[] F;
    return 0;
}