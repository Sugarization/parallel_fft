#ifndef _TEST_FFT_H__
#define _TEST_FFT_H__

#include "fft_base.h"
#include "fft_parallel.h"
#include "timer.h"

#include <omp.h>
#include <cstdlib>
#include <cstdio>

enum class FFT_Type{naive, iter, cooley, embed};

double fft_test(NComplex *F, NComplex *x, indexT N, FFT_Type type, int nThreads, int inv = -1);
void sanityCheck(FFT_Type t1, FFT_Type t2, indexT N, int inv);
void Speedtest(FFT_Type type, indexT N, int rep, int nThreads);

void NMMSE(NComplex *y, NComplex *x, indexT N, double &max, double &mean);

void copyComplex(NComplex *y, NComplex *x, indexT N);
void printArray(NComplex *x, indexT N);

class Generator
{
    double rand01();

public:
    Generator();
    void linear(NComplex *x, indexT N, NComplex step);
    void sine(NComplex *x, indexT N, double omega);
    void random(NComplex *x, indexT N);
};

#endif