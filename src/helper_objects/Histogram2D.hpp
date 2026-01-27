#pragma once

#include <string>

class Histogram2D
{
    private:

    public:
        int numBinsX;
        double minValX;
        double maxValX;
        double binSizeX;
        double* binWallsX;

        int numBinsY;
        double minValY;
        double maxValY;
        double binSizeY;
        double* binWallsY;
        
        // Both Y major 2D arrays
        int* counts_par;
        int* counts_perp;

        int outOfBoundsCount = 0;

        Histogram2D(int inNumBinsX, int inNumBinsY, double inMinX, double inMaxX, double inMinY, double inMaxY, 
            bool logBinsX=false, bool logBinsY=false);

        ~Histogram2D();

        void addVal(double xValue, double yValue, int pol);

        void combineData(Histogram2D& otherHist);

        void overrideVal(double xValue, double yValue, int pol, int val);

        void exportToFile(const std::string& name, const std::string& folder = "None");
};