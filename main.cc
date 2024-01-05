#include "fft_base.h"
#include <cstdlib>
#include <ctime>
#include <cstdio>

const indexT N = 1 << 15;
const T Fs = 64;
NComplex x[N], F[N];

int main()
{
    srand(time(nullptr));
    T f = 4;
    for (indexT i = 0; i < N; ++ i) {
        x[i].re = sin(2 * PI * f / Fs * i);
        x[i].re += i;
        x[i].im = 0;
    }
    naiveDFT(F, x, N);
    iterativeFFT(x, N);
    for (indexT i = 0; i < N; ++ i) {
        if ((x[i] - F[i]).abs() > 1e-6) {
            puts("check false");
        }
    }
    
    return 0;
}