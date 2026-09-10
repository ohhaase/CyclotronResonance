#include "AvgOutput.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include "global_vars.hpp"

AvgOutput::AvgOutput(int inNumBinsX, int inNumBinsY, double inMinX, double inMaxX, double inMinY, double inMaxY, 
    bool logBinsX, bool logBinsY)
{
    // Stick input values into respective variables
    numBinsX = inNumBinsX;
    minValX = inMinX;
    maxValX = inMaxX;

    numBinsY = inNumBinsY;
    minValY = inMinY;
    maxValY = inMaxY;

    // Initialize the X bin walls
    binWallsX.resize(numBinsX + 1);
    binWallsX[0] = minValX;

    if (logBinsX)
    {
        double logBinSizeX = (log10(maxValX) - log10(minValX)) / numBinsX;

        for (int i = 0; i < numBinsX; i++)
        {
            binWallsX[i+1] = pow(10, log10(binWallsX[i]) + logBinSizeX);
        }
    }
    else
    {
        double binSizeX = (maxValX - minValX) / numBinsX;

        for (int i = 0; i < numBinsX; i++)
        {
            binWallsX[i+1] = binWallsX[i] + binSizeX;
        }
    }

    // Initialize the Y bin walls
    binWallsY.resize(numBinsY + 1);
    binWallsY[0] = minValY;

    if (logBinsY)
    {
        double logBinSizeY = (log10(maxValY) - log10(minValY)) / numBinsY;

        for (int i = 0; i < numBinsY; i++)
        {
            binWallsY[i+1] = pow(10, log10(binWallsY[i]) + logBinSizeY);
        }
    }
    else
    {
        double binSizeY = (maxValY - minValY) / numBinsY;

        for (int i = 0; i < numBinsY; i++)
        {
            binWallsY[i+1] = binWallsY[i] + binSizeY;
        }
    }

    // Initialize the count arrays
    counts.resize(numBinsX * numBinsY);
    total_nrg.resize(numBinsX * numBinsY);
    total_theta.resize(numBinsX * numBinsY);
    total_escape_count.resize(numBinsX * numBinsY);
    total_polarization.resize(numBinsX * numBinsY);

}


void AvgOutput::addVal(PhotonState initPhoton, PhotonState photon)
{
    // Loop over all bins to see which it goes in
    // If a value doesn't fall within these ranges, then keep track of how many

    // For avg data, x is always nrg and y is always theta

    for (int i = 0; i < numBinsX; i++)
    {
        if ((binWallsX[i] <= initPhoton.omega) && (initPhoton.omega < binWallsX[i+1]))
        {
            for (int j = 0; j < numBinsY; j++)
            {
                if ((binWallsY[j] <= initPhoton.theta) && (initPhoton.theta < binWallsY[j+1]))
                {
                    counts[i*numBinsY + j] += 1;
                    total_nrg[i*numBinsY + j] += photon.omega;
                    total_theta[i*numBinsY + j] += photon.theta;
                    total_escape_count[i*numBinsY + j] += static_cast<double>(photon.numScatterings);
                    total_polarization[i*numBinsY + j] += static_cast<double>(photon.polarization);
                    return; // Stop if we found it
                }
            }
        }
    }
    
    // If we get to the end and haven't found it, add one to the ticker
    outOfBoundsCount += 1;
}


void AvgOutput::exportToFile(const std::string& name, const std::string& folder)
{
    std::string fileName = name + ".csv";

    if (folder != "None")
    {
        // Note: folder must already exist!
        fileName = folder + "/" + fileName;
    }

    std::ofstream file(fileName);

    /* 
    File has the following structure:
    | omega | theta | avg_nrg | avg_theta | avg_scatterings | avg_polarization | 0
    0 is there to make files legacy compatible
    */

    for (int i = 0; i < numBinsX; i++)
    {
        for (int j = 0; j < numBinsY; j++)
        {
            int ind = i*numBinsY + j;

            file << 0.5*(binWallsX[ind] + binWallsX[ind + numBinsY]) << "," << 
                0.5*(binWallsY[ind] + binWallsY[ind + 1]) << "," << 
                total_nrg[ind] / counts[ind] << "," <<
                total_theta[ind] / counts[ind] << "," <<
                total_escape_count[ind] / counts[ind] << "," <<
                total_polarization[ind] / counts[ind] << "," <<
                "0\n";
        }
    }

    file.close();
}