#ifndef _TIMER_H__
#define _TIMER_H__

#include <ctime>
#include <cstring>
#include <cstdint>
#include <cstdio>

class Timer
{
    timespec t0, t1;
    char *name;
public:
    Timer(const char *);
    ~Timer();
    void tick();
    void tock(bool verbose = true);
    double difftime();
};

#endif