#include "fft_base.h"
#ifndef NDEBUG
#include <cstdio>
#endif

NComplex::NComplex(): re(.0), im(.0) {}
NComplex::NComplex(T r, T i): re(r), im(i) {}

// NComplex NComplex::operator + (const NComplex& x) const {
//     __m128d a, b;
//     a = _mm_load_pd(&re);
//     b = _mm_load_pd(&x.re);
//     a = _mm_add_pd(a, b);
//     return *(NComplex*)&a;
// }
// NComplex NComplex::operator - (const NComplex& x) const {
//     __m128d a, b;
//     a = _mm_load_pd(&re);
//     b = _mm_load_pd(&x.re);
//     a = _mm_sub_pd(a, b);
//     return *(NComplex*)&a;
// }
NComplex NComplex::operator * (const NComplex& x) const {
    __m128d s, o, sr, si, oswap, zero;
    s = _mm_load_pd(&re);
    o = _mm_load_pd(&x.re);
    sr = _mm_shuffle_pd(s, s, 0);                   // self: real, real
    si = _mm_shuffle_pd(s, s, 3);                   // self: imag, imag
    oswap = _mm_shuffle_pd(o, o, 1);                // other: imag, real
    zero = _mm_set1_pd(0.0);
    oswap = _mm_addsub_pd(zero, oswap); // other: -imag, real
    o = _mm_mul_pd(o, sr);
    oswap = _mm_mul_pd(oswap, si);
    o = _mm_add_pd(o, oswap);
    return *(NComplex*)&o;
}
NComplex NComplex::operator + (const NComplex& x) const { return {re + x.re, im + x.im}; }
NComplex NComplex::operator - (const NComplex& x) const { return {re - x.re, im - x.im}; }
// NComplex NComplex::operator * (const NComplex& x) const { return {re * x.re - im * x.im, re * x.im + im * x.re}; }
NComplex NComplex::operator * (T k) const { return {re * k, im * k}; }
NComplex NComplex::operator / (const NComplex& x) const { T nsq = normSq(); return {re / nsq, -im / nsq}; }
T NComplex::normSq() const { return re * re + im * im; }
T NComplex::abs() const { return sqrt(normSq()); }
T NComplex::phase() const { return atan2(im, re); }
NComplex expJ(T theta) { return {cos(theta), sin(theta)}; }
NComplex NComplex::j() const { return {-im, re}; }

void zeroC(NComplex *x, indexT N)
{
    memset(x, 0, sizeof(NComplex) * N);
}

indexT bitReverse(indexT x_in, int K, int grain)
{
    static const uint64_t rMask[8] = {0x5555555555555555ull, 0x3333333333333333ull, 0x0f0f0f0f0f0f0f0full, 0x00ff00ff00ff00ffull, 0x0000ffff0000ffffull};
    static const uint64_t lMask[8] = {0xaaaaaaaaaaaaaaaaull, 0xccccccccccccccccull, 0xf0f0f0f0f0f0f0f0ull, 0xff00ff00ff00ff00ull, 0xffff0000ffff0000ull};
    uint64_t x = x_in;
    // if (grain == 0) 
    x = ((x & rMask[0]) << 1) | ((x & lMask[0]) >> 1);
    x = ((x & rMask[1]) << 2) | ((x & lMask[1]) >> 2);
    x = ((x & rMask[2]) << 4) | ((x & lMask[2]) >> 4);
    x = ((x & rMask[3]) << 8) | ((x & lMask[3]) >> 8);
    x = ((x & rMask[4]) << 16) | ((x & lMask[4]) >> 16);
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

void bitReverseSwap(NComplex *x, indexT N)
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
}

void butterfly(NComplex *x, indexT N, indexT m, int inv)
{
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

void columnButterfly(NComplex *x, indexT j, indexT N1, indexT N2, indexT m, int inv)
{
    int K2 = MSB(N2);
    NComplex w1 = expJ(inv * PI / m);

    for (indexT s = 0; s < N1; s += 2 * m) {
        
        NComplex w = {1, 0};
        for (indexT i = 0; i < m; ++ i) {
            indexT i1 = P2RC(s + i, j, K2);
            indexT i2 = P2RC(s + i + m, j, K2);
            
            NComplex u = x[i1], v = x[i2] * w;
            x[i1] = u + v;
            x[i2] = u - v;
            w = w * w1;
        }
    }
}

void iterativeFFT(NComplex *x, indexT N, int inv)
{
    bitReverseSwap(x, N);
    for (indexT m = 1; m < N; m <<= 1) {
        butterfly(x, N, m, inv);
    }
}

void iterativeColumnFFT(NComplex *x, indexT cl, indexT N1, indexT N2, int inv)
{
    int K1 = MSB(N1), K2 = MSB(N2);
    for (indexT i = 0; i < N1; i ++) {
        indexT j = bitReverse(i, K1);
        if (i < j) {
            NComplex tempor = x[P2RC(i, cl, K2)];
            x[P2RC(i, cl, K2)] = x[P2RC(j, cl, K2)];
            x[P2RC(j, cl, K2)] = tempor;
        }
    }
    for (indexT m = 1; m < N1; m <<= 1) {
        columnButterfly(x, cl, N1, N2, m, inv);
    }
}

void iterativeRowFFT(NComplex *x, indexT ro, indexT N1, indexT N2, int inv)
{
    iterativeFFT(x + ro * N2, N2, inv);
}