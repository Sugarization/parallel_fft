#include "timer.h"

void Timer::tick()
{
    clock_gettime(CLOCK_MONOTONIC, &t0);
}

double Timer::difftime() {
    uint64_t diff = (t1.tv_sec - t0.tv_sec) * 1000000000ull + (t1.tv_nsec - t0.tv_nsec);
    return diff / 1e9;
}

void Timer::tock(const char *name)
{
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (name) {
        printf("Procedure %s elapsed %.7f s.\n", name, difftime());
    }
}
