#include "Timer.hpp"

#include <chrono>

Timer::Timer()
{
    startTime = std::chrono::high_resolution_clock::now();
}


Timer::~Timer()
{

}


double Timer::elapsedNano()
{
    endTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::nano> dT = endTime - startTime;

    return dT.count();
}


double Timer::elapsedMicro()
{
    endTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::micro> dT = endTime - startTime;

    return dT.count();
}


double Timer::elapsedMilli()
{
    endTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> dT = endTime - startTime;

    return dT.count();
}


double Timer::elapsedSec()
{
    endTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dT = endTime - startTime;

    return dT.count();
}