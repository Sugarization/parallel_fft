#include "tests.h"

const char *name(FFT_Type type)
{
    switch (type)
    {
    case FFT_Type::iter:
        return "Iterative Cooley-Tucky";
        
    case FFT_Type::cooley:
        return "Embarrasingly Parallel Cooley-Tucky";

    case FFT_Type::embed:
        return "Embedded into 2D";

    case FFT_Type::naive:
        return "Naive DFT";
    }
    return "";
}

void copyComplex(NComplex *y, NComplex *x, indexT N)
{
    memcpy(y, x, sizeof(NComplex) * N);
}

void printArray(NComplex *x, indexT N)
{
    for (indexT i = 0; i < N; ++ i) {
        printf("(%.3f, %.3f) ", x[i].re, x[i].im);
    }
    puts("");
}

void sanityCheck(FFT_Type t1, FFT_Type t2, indexT N, int inv)
{
    auto F = new NComplex[N];
    auto F0 = new NComplex[N];
    auto x = new NComplex[N];
    auto x0 = new NComplex[N];

    Generator g;
    g.random(x, N);
    copyComplex(x0, x, N);

    fft_test(F, x, N, t1, -1, inv);
    fft_test(F0, x0, N, t2, -1, inv);
    double max, mean;
    NMMSE(F0, F, N, max, mean);
    printf("%.12f %.12f\n", max, mean);
    // printArray(F, N);
    // printArray(F0, N);
    delete[] F;
    delete[] x;
    delete[] F0;
    delete[] x0;
}

void Speedtest(FFT_Type type, indexT N, int rep, int nThreads)
{
    auto F = new NComplex[N];
    auto x = new NComplex[N];

    Generator g;
    g.random(x, N);

    double total = .0;
    for (int r = 0; r < rep; ++ r) {
        total += fft_test(F, x, N, type, nThreads);
    }

    printf("Method \'%s\' averages %.7f s on N = %lu with max %d cores (repeated %d times).\n", name(type), total / rep, N, nThreads, rep);

    delete[] F;
    delete[] x;
}

double fft_test(NComplex *F, NComplex *x, indexT N, FFT_Type type, int nThreads, int inv)
{
    // Test a fft function of _type_ using input _x_ and outputs into _F_
    // Warning: this function possibly modifies _x_.
    // returns running time
    if (nThreads > 0) {
        omp_set_num_threads(nThreads);
    }

    Timer tm;
    tm.tick();
   
    switch (type)
    {
        case FFT_Type::iter:
        copyComplex(F, x, N);
        iterativeFFT(F, N, inv);
        break;

        case FFT_Type::cooley:
        directParallelFFT(F, x, N, inv);
        break;

        case FFT_Type::embed:
        embeddingParallelColumnCopyFFT(F, x, N, inv);
        break;

        case FFT_Type::naive:
        naiveDFT(F, x, N, inv);
        break;

    default:
        puts("FFT Type error");
        break;
    }

    tm.tock();
    return tm.difftime();
}

void NMMSE(NComplex *y, NComplex *x, indexT N, double &max, double &mean)
{
    // normalized mean and max squared error (for fft outputs)
    max = 0.;
    mean = 0.;
    for (indexT i = 0; i < N; ++ i) {
        double error = (x[i] - y[i]).normSq() / N;
        if (error > max) max = error;
        mean += error;
    }
    mean /= N;
}

Generator::Generator()
{
    srand(time(nullptr));
}

void Generator::linear(NComplex *x, indexT N, NComplex step)
{
    x[0] = {0, 0};
    for (indexT i = 0; i < N; ++ i) {
        x[i] = x[i - 1] + step;
    }
}

void Generator::sine(NComplex *x, indexT N, double omega)
{
    for (indexT i = 0; i < N; ++ i) {
        x[i] = {sin(omega * i), 0};
    }
}

void Generator::random(NComplex *x, indexT N)
{
    for (indexT i = 0; i < N; ++ i) {
        x[i] = {rand01() - 0.5, rand01() - 0.5};
    }
}

double Generator::rand01()
{
    return rand() / (double)RAND_MAX;
}