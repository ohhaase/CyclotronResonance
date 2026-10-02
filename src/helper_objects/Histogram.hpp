#pragma once

#include <string>
#include <cstdint>
#include <vector>

class Histogram
{
    private:


    public:
        int numBins;
        double minVal;
        double maxVal;
        std::vector<double> binWalls;
        std::vector<uint64_t> counts_par;
        std::vector<uint64_t> counts_perp;

        int outOfBoundsCount = 0;

        Histogram(int inNumBins, double inMin, double inMax, bool logBins=false);

        void addVal(double value, int pol);

        void exportToFile(const std::string& name, const std::string& folder = "None");

        void combineData(Histogram& otherHist);
};