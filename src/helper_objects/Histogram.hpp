#pragma once

#include <string>

class Histogram
{
    private:


    public:
        int numBins;
        double minVal;
        double maxVal;
        double* binWalls;
        int* counts_par;
        int* counts_perp;

        int outOfBoundsCount = 0;

        Histogram(int inNumBins, double inMin, double inMax, bool logBins=false);

        ~Histogram();

        void addVal(double value, int pol);

        void exportToFile(const std::string& name, const std::string& folder = "None");

        void combineData(Histogram& otherHist);
};