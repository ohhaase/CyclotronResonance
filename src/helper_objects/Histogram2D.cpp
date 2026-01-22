#include "Histogram2D.hpp"

#include <cmath>
#include <string>
#include <fstream>

Histogram2D::Histogram2D(int inNumBinsX, int inNumBinsY, double inMinX, double inMaxX, double inMinY, double inMaxY, 
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
    binWallsX = new double[numBinsX + 1]{};
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
    binWallsY = new double[numBinsY + 1]{};
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
    counts_par = new int[numBinsX * numBinsY]{}; // Initialized to zero
    counts_perp = new int[numBinsX * numBinsY]{}; // Initialized to zero

}


Histogram2D::~Histogram2D()
{
    delete[] binWallsX;
    delete[] binWallsY;
    delete[] counts_par;
    delete[] counts_perp;
}


void Histogram2D::addVal(double xValue, double yValue, int pol)
{
    // Figure out which count we add to based on polarization
    int* thisCount;

    switch (pol)
    {
        case -1: // Average polarization
        {
            // In this case, we just treat parallel count as the total count
            thisCount = counts_par;
            break;
        }
        case 0: // Parallel
        {
            thisCount = counts_par;
            break;
        }
        case 1: // Perp
        {
            thisCount = counts_perp;
            break;
        }

    }

    // Loop over all bins to see which it goes in
    // If a value doesn't fall within these ranges, then keep track of how many
    for (int i = 0; i < numBinsX; i++)
    {
        if ((binWallsX[i] <= xValue) && (xValue < binWallsX[i+1]))
        {
            for (int j = 0; j < numBinsY; j++)
            {
                if ((binWallsY[j] <= yValue) && (yValue < binWallsY[j+1]))
                {
                    thisCount[i*numBinsY + j] += 1;
                    return; // Stop if we found it
                }
            }
        }
    }
    
    // If we get to the end and haven't found it, add one to the ticker
    outOfBoundsCount += 1;
}


void Histogram2D::overrideVal(double xValue, double yValue, int pol, int val)
{
    int* thisCount = counts_par;

    switch (pol)
    {
        case -1: // Average polarization
        {
            // In this case, we just treat parallel count as the total count
            thisCount = counts_par;
            break;
        }
        case 0: // Parallel
        {
            thisCount = counts_par;
            break;
        }
        case 1: // Perp
        {
            thisCount = counts_perp;
            break;
        }

    }

    // Loop over all bins to see which it goes in
    // Since this will usually be used as a 2D array, these should be unique if cell centered
    for (int i = 0; i < numBinsX; i++)
    {
        if ((binWallsX[i] <= xValue) && (xValue < binWallsX[i+1]))
        {
            for (int j = 0; j < numBinsY; j++)
            {
                if ((binWallsY[j] <= yValue) && (yValue < binWallsY[j+1]))
                {
                    thisCount[i*numBinsY + j] = val;
                    return; // Stop if we found it
                }
            }
        }
    }
}


void Histogram2D::exportToFile(const std::string& name, const std::string& folder)
{
    std::string fileName = "hist2D_" + name + ".csv";

    if (folder != "None")
    {
        // Note: folder must already exist!
        fileName = folder + "/" + fileName;
    }

    std::ofstream file(fileName);

    /* 
    File has the following structure:
    | -- | x0 | -- | x1 | -- | .... | xn | -- |xn+1| -- |
    | y0 | pl | pr | pl | pr | .... | pr | pl | -- | -- |
    | y1 | pl | pr | pl | pr | .... | pr | pl | -- | -- |
        .                         .                    .
        .                         .                    .
    | yn | pl | pr | pl | pr | .... | pr | pl | -- | -- |
    |yn+1| -- | -- | -- | -- | .... | -- | -- | -- | -- |
    
    Extra zeroes (--) are necessary for easy python parsing
    */

    // Write a zero as a spacer first
    file << "0,";

    // Write x wall values to file first
    for (int i = 0; i < numBinsX + 1; i++)
    {
        // Need to include zeros every other spot to account for both counts
        file << binWallsX[i] << ",0,";
    }

    // Newline
    file << "\n";

    // Write y walls and count data in a big loop
    for (int j = 0; j < numBinsY; j++)
    {
        // Write y wall value
        file << binWallsY[j] << ",";

        for (int i = 0; i < numBinsX; i++)
        {
            file << counts_par[i*numBinsY + j] << "," << counts_perp[i*numBinsY + j] << ",";
        }

        // Spacers at the end of a line
        file << "0,0,\n";
    }
    
    // Add last y value
    file << binWallsY[numBinsY] << ",";

    // Add zero padding at bottom
    for (int i = 0; i < (numBinsX+1)*2; i++)
    {
        file << "0,";
    }

    file.close();
}