#include "transpose.h"

void parallelNaiveTranspose(NComplex *y, NComplex *x, indexT N1, indexT N2)
{
    int K1 = MSB(N1), K2 = MSB(N2);
    #pragma omp parallel for 
    for (indexT i = 0; i < N1; ++ i){
        for (indexT j = 0; j < N2; ++ j) {
            y[P2RC(j, i, K1)] = x[P2RC(i, j, K2)];
        }
    }
}

inline void transposeBlock(NComplex *y, NComplex *x, indexT K1, indexT K2)
{
    for (indexT i = 0; i < BS; ++ i) {
        for (indexT j = 0; j < BS; ++ j) {
            y[P2RC(j, i, K1)] = x[P2RC(i, j, K2)];
        }
    }
}

void parallelTranspose(NComplex *y, NComplex *x, indexT N1, indexT N2)
{
// This function assumes that N1 and N2 are powers of 2.
    if (N1 <= 4 * BS && N2 <= 4 * BS) {
        parallelNaiveTranspose(y, x, N1, N2);
        return;
    }

    int K1 = MSB(N1), K2 = MSB(N2);
    #pragma omp parallel for 
    for (indexT i = 0; i < N1; i += BS) {
        for (indexT j = 0; j < N2; j += BS) {
            transposeBlock(&y[P2RC(j, i, K1)], &x[P2RC(i, j, K2)], K1, K2);
        }
    }
}