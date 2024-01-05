#ifndef _FFT_PARALLEL_H_
#define _FFT_PARALLEL_H_

#include "fft_base.h"
#include <omp.h>

void directParallelFFT(NComplex *F, NComplex *x, indexT N, int inv = -1);

#endif 