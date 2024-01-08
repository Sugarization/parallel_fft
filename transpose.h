#ifndef _TRANSPOSE_H_
#define _TRANSPOSE_H_

#include <omp.h>
#include <immintrin.h>
#include "fft_base.h"

#define BS 16 // block size

void parallelNaiveTranspose(NComplex *y, NComplex *x, indexT N1, indexT N2);
void parallelTranspose(NComplex *y, NComplex *x, indexT N1, indexT N2);

#endif