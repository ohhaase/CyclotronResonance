#pragma once

#include <vector>
#include <string>

class NumericalDistb
{
    // Uses two arrays, xVals and probVals, to make a samplable distribution from discrete values
    private:
        int N;
        std::vector<double> xVals;
        std::vector<double> probVals;

        bool logSpacing;

        // Needs to support log spacing.
        // If creating with logSpacing = true, these correspond to x values in log space
        double xMin;
        double xMax;
        double dx;

    public:
        NumericalDistb(const std::string& filePath, bool inLogSpacing);

        ~NumericalDistb();

        double evaluateDistribution(double x);
};