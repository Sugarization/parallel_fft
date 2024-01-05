#include "fft_base.h"
#ifndef NDEBUG
#include <cstdio>
#endif

NComplex::NComplex(): re(.0), im(.0) {}
NComplex::NComplex(T r, T i): re(r), im(i) {}

NComplex NComplex::operator + (const NComplex& x) const { return {re + x.re, im + x.im}; }
NComplex NComplex::operator - (const NComplex& x) const { return {re - x.re, im - x.im}; }
NComplex NComplex::operator * (const NComplex& x) const { return {re * x.re - im * x.im, re * x.im + im * x.re}; }
NComplex NComplex::operator * (T k) const { return {re * k, im * k}; }
NComplex NComplex::operator / (const NComplex& x) const { T nsq = normSq(); return {re / nsq, -im / nsq}; }
T NComplex::normSq() const { return re * re + im * im; }
T NComplex::abs() const { return sqrt(normSq()); }
T NComplex::phase() const { return atan2(im, re); }
NComplex expJ(T theta) { return {cos(theta), sin(theta)}; }

void zeroC(NComplex *x, indexT N)
{
    memset(x, 0, sizeof(NComplex) * N);
}

indexT bitReverse(indexT x_in, int K)
{
    static const uint64_t rMask[8] = {0x5555555555555555ull, 0x3333333333333333ull, 0x0f0f0f0f0f0f0f0full, 0x00ff00ff00ff00ffull, 0x0000ffff0000ffffull};
    static const uint64_t lMask[8] = {0xaaaaaaaaaaaaaaaaull, 0xccccccccccccccccull, 0xf0f0f0f0f0f0f0f0ull, 0xff00ff00ff00ff00ull, 0xffff0000ffff0000ull};
    uint64_t x = x_in;
    x = ((x & rMask[0]) << 1) | ((x & lMask[0]) >> 1);
    x = ((x & rMask[1]) << 2) | ((x & lMask[1]) >> 2);
    x = ((x & rMask[2]) << 4) | ((x & lMask[2]) >> 4);
    x = ((x & rMask[3]) << 8) | ((x & lMask[3]) >> 8);
    if (K <= 16) return x >> (16 - K);

    x = ((x & rMask[4]) << 8) | ((x & lMask[4]) >> 8);
    x = ((x & rMask[5]) << 16) | ((x & lMask[5]) >> 16);
    if (K <= 32) return x >> (32 - K);

    x = (x << 32 | x >> 32);
    return x >> (64 - K);
}

void bitReverse(indexT *x, indexT N)
{
    int K = MSB(N);
    for (int i = 0; i < N; ++ i) {
        x[i] = bitReverse(x[i], K);
    }
}

void naiveDFT(NComplex *F, NComplex *x, indexT N, int inv)
{
    for (indexT k = 0; k < N; ++ k) {
        F[k] = {0, 0};
        NComplex wNk = expJ(inv * 2 * PI * k / N);
        NComplex w = {1, 0};
        for (indexT n = 0; n < N; ++ n) {
            F[k] = F[k] + x[n] * w;
            w = w * wNk;
        }
    }
}

void iterativeFFT(NComplex *x, indexT N, int inv)
{
    int K = MSB(N);
    for (indexT i = 0; i < N; i ++) {
        indexT j = bitReverse(i, K);
        if (i < j) {
            NComplex tempor = x[i];
            x[i] = x[j];
            x[j] = tempor;
        }
    }
    for (indexT m = 1; m < N; m <<= 1) {
        NComplex w1 = expJ(inv * PI / m);

        for (indexT s = 0; s < N; s += 2 * m) {
            
            NComplex w = {1, 0};
            for (indexT i = 0; i < m; ++ i) {
                NComplex u = x[s + i];
                NComplex v = x[s + m + i] * w;
                
                x[s + i] =     u + v;
                x[s + i + m] = u - v;
                w = w * w1;
            }
        }
    }
}