#ifndef _TIMER_H__
#define _TIMER_H__

#include <ctime>
#include <cstring>
#include <cstdint>
#include <cstdio>

class Timer
{
    timespec t0, t1;
public:
    void tick();
    void tock(const char* name = nullptr);
    double difftime();
};

#endif