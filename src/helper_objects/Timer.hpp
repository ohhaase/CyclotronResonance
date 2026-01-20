#pragma once

#include <chrono>

class Timer
{
    private:
        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;
    
    public:
        Timer();

        ~Timer();

        double elapsedNano();

        double elapsedMicro();

        double elapsedMilli();

        double elapsedSec();
};