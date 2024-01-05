#ifndef _FFT_BASE_H_
#define _FFT_BASE_H_

#include <cmath>
#include <cstdint>
#include <cstring>

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
    T normSq() const;
    T abs() const;
    T phase() const;
};

NComplex expJ(T theta);

#define PI M_PI
#define MSB(n) (sizeof(indexT) * 8 - 1 - __builtin_clzll((unsigned long long)n))

indexT bitReverse(indexT x, int K);
void   bitReverse(indexT *x, indexT N);

void zeroC(NComplex *x, indexT N);
void naiveDFT(NComplex *F, NComplex *x, indexT N, int inv = -1);
void iterativeFFT(NComplex *x, indexT N, int inv = -1);

#endif