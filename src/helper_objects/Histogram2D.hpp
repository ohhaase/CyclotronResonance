#pragma once

#include <string>
#include <vector>

class Histogram2D
{
    private:

    public:
        int numBinsX;
        double minValX;
        double maxValX;
        double binSizeX;
        std::vector<double> binWallsX;

        int numBinsY;
        double minValY;
        double maxValY;
        double binSizeY;
        std::vector<double> binWallsY;
        
        // Both Y major 2D arrays
        std::vector<int> counts_par;
        std::vector<int> counts_perp;

        int outOfBoundsCount = 0;

        Histogram2D(int inNumBinsX, int inNumBinsY, double inMinX, double inMaxX, double inMinY, double inMaxY, 
            bool logBinsX=false, bool logBinsY=false);

        void addVal(double xValue, double yValue, int pol);

        void combineData(Histogram2D& otherHist);

        void overrideVal(double xValue, double yValue, int pol, int val);

        void exportToFile(const std::string& name, const std::string& folder = "None");
};