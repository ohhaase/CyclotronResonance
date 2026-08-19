#pragma once

#include <string>
#include <cstdint>

class Histogram
{
    private:


    public:
        int numBins;
        double minVal;
        double maxVal;
        double* binWalls;
        uint64_t* counts_par;
        uint64_t* counts_perp;

        int outOfBoundsCount = 0;

        Histogram(int inNumBins, double inMin, double inMax, bool logBins=false);

        ~Histogram();

        void addVal(double value, int pol);

        void exportToFile(const std::string& name, const std::string& folder = "None");

        void combineData(Histogram& otherHist);
};