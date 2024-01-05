#include "fft_base.h"
#include "fft_parallel.h"
#include <cstdlib>
#include <ctime>
#include <cstdio>

const indexT N = 1 << 25;
const T Fs = 64;
NComplex x[N], x1[N], F[N];

int main()
{
    srand(time(nullptr));
    T f = 4;
    for (indexT i = 0; i < N; ++ i) {
        x[i].re = sin(2 * PI * f / Fs * i);
        x[i].re += i;
        x[i].im = 0;
        x1[i] = x[i];
    }
    // naiveDFT(F, x, N);
    directParallelFFT(F, x, N);
    #ifndef NDEBUG
    iterativeFFT(x1, N);
    for (indexT i = 0; i < N; ++ i) {
        T error = (x1[i] - F[i]).abs();
        if (error > 1e-6) {
            printf("Error check failed: %.8f\n", error);
        }
    }
    #endif
    
    return 0;
}