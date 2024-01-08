#ifndef _FFT_BASE_H_
#define _FFT_BASE_H_

#include <cmath>
#include <cstdint>
#include <cstring>
#include <immintrin.h>

#define P2RC(i, j, K) (((i) << (K)) | (j))

typedef double T;
typedef uint64_t indexT;

struct NComplex
{
    T re;
    T im;
    NComplex();
    NComplex(T, T);
    NComplex operator + (const NComplex&) const;
    NComplex operator - (const NComplex&) const;
    NComplex operator * (const NComplex&) const;
    NComplex operator * (T) const;
    NComplex operator / (const NComplex&) const;
    NComplex j() const;
    T normSq() const;
    T abs() const;
    T phase() const;
};

NComplex expJ(T theta);

#define PI M_PI
#define MSB(n) (sizeof(indexT) * 8 - 1 - __builtin_clzll((unsigned long long)n))

indexT bitReverse(indexT x, int K, int grain = 0);
void   bitReverse(indexT *x, indexT N);

void zeroC(NComplex *x, indexT N);
void naiveDFT(NComplex *F, NComplex *x, indexT N, int inv = -1);
void iterativeFFT(NComplex *x, indexT N, int inv = -1);

void iterativeRowFFT(NComplex *x, indexT ro, indexT N1, indexT N2, int inv);
void iterativeColumnFFT(NComplex *x, indexT cl, indexT N1, indexT N2, int inv);

#endif