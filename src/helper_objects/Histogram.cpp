#include "Histogram.hpp"

#include <cmath>
#include <string>
#include <fstream>
#include <iostream>


Histogram::Histogram(int inNumBins, double inMin, double inMax, bool logBins)
{
    // Stick input values into respective variables
    numBins = inNumBins;
    minVal = inMin;
    maxVal = inMax;

    // Initialize the bin walls
    binWalls = new double[numBins + 1]{};
    binWalls[0] = minVal;

    if (logBins)
    {
        double logBinSize = (log10(maxVal) - log10(minVal)) / numBins;

        for (int i = 0; i < numBins; i++)
        {
            binWalls[i+1] = pow(10, log10(binWalls[i]) + logBinSize);
        }
    }
    else
    {
        double binSize = (maxVal - minVal) / numBins;

        for (int i = 0; i < numBins; i++)
        {
            binWalls[i+1] = binWalls[i] + binSize;
        }
    }

    // Initialize the other arrays
    counts_par = new int[numBins]{}; // Initialized to zero
    counts_perp = new int[numBins]{}; // Initialized to zero
}


Histogram::~Histogram()
{
    delete[] binWalls;
    delete[] counts_par;
    delete[] counts_perp;
}


void Histogram::addVal(double value, int pol)
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
    // Could speed this up with a binary search? Hardly feels necessary lol
    // If a value doesn't fall within these ranges, then... oh well? keep track of how many?
    for (int i = 0; i < numBins; i++)
    {
        if ((binWalls[i] <= value) && (value < binWalls[i+1]))
        {
            thisCount[i] += 1;
            return; // Stop if we found it
        }
    }
    
    // If we get to the end and haven't found it, add one to the ticker
    outOfBoundsCount += 1;
}


void Histogram::exportToFile(const std::string& name, const std::string& folder)
{
    std::string fileName = "hist_" + name + ".csv";

    if (folder != "None")
    {
        // Note: folder must already exist!
        fileName = folder + "/" + fileName;
    }

    std::ofstream file(fileName);

    // Write wall values to file
    for (int i = 0; i < numBins + 1; i++)
    {
        file << binWalls[i] << ",";
    }

    // Newline
    file << "\n";
    
    // Write parallel counts to file
    for (int i = 0; i < numBins; i++)
    {
        file << counts_par[i] << ",";
    }

    // Need to pad with a zero, because binWalls is one bigger
    file << "0,\n";

    // Write perp counts to file
    for (int i = 0; i < numBins; i++)
    {
        file << counts_perp[i] << ",";
    }

    // Need to pad with a zero, because binWalls is one bigger
    file << "0,";

    file.close();
}


void Histogram::combineData(Histogram& otherHist)
{
    // Combine count data from other histogram
    // Assumes all values are the same!! We don't stop if they aren't but we do flag it
    if ((numBins != otherHist.numBins) || (minVal != otherHist.minVal) || (maxVal != otherHist.maxVal))
    {
        std::cout << "Incompatible histograms! minVal: " << minVal << " maxVal: " << maxVal;
    }

    // Simply loop over count arrays and add them together
    for (int i = 0; i < numBins; i++)
    {
        counts_par[i] += otherHist.counts_par[i];
        counts_perp[i] += otherHist.counts_perp[i];
    }
}