#ifndef _TEST_FFT_H__
#define _TEST_FFT_H__

#include "fft_base.h"
#include "fft_parallel.h"
#include "timer.h"

#include <omp.h>
#include <cstdlib>
#include <cstdio>

enum class FFT_Type{iter, ct, embed};

void fft_test(NComplex *F, NComplex *x, indexT N, FFT_Type type, int nThreads, int inv = -1);

#endif