#include "fft_base.h"
#include "fft_parallel.h"
#include "timer.h"

#include <omp.h>
#include <cstdlib>
#include <cstdio>

const indexT N = 1 << 25;
const T Fs = 64;


int main()
{
    omp_set_num_threads(16);
    NComplex *x = new NComplex[N];
    NComplex *F = new NComplex[N];
    
    srand(time(nullptr));
    T f = 10;
    for (indexT i = 0; i < N; ++ i) {
        // x[i].re = sin(2 * PI * f / Fs * i);
        x[i].re = i;
        x[i].im = 0;
    }
    #ifndef FFT_CHECK
    NComplex *x1 = new NComplex[N];
    for (indexT i = 0; i < N; ++ i) {
        x1[i] = x[i];
        // printf("%.4f %.4f\n", x1[i].re, x1[i].im);
    }
    #endif
    Timer timer("Direct Parallel");
    timer.tick();
    directParallelFFT(F, x, N);
    timer.tock();

    #ifdef FFT_CHECK
    
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