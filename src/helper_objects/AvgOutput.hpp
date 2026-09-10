#pragma once

#include <string>
#include <vector>

#include "global_vars.hpp"

class AvgOutput
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
        
        // Y major 2D arrays
        std::vector<double> total_nrg;
        std::vector<double> total_theta;
        std::vector<double> total_escape_count;
        std::vector<double> total_polarization;
        std::vector<int> counts;

        int outOfBoundsCount = 0;

        AvgOutput(int inNumBinsX, int inNumBinsY, double inMinX, double inMaxX, double inMinY, double inMaxY, 
            bool logBinsX=false, bool logBinsY=false);

        void addVal(PhotonState initPhoton, PhotonState photon);

        void exportToFile(const std::string& name, const std::string& folder = "None");
};