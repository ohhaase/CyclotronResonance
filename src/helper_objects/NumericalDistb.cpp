#include "NumericalDistb.hpp"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>

NumericalDistb::NumericalDistb(const std::string& filePath, bool inLogSpacing)
{
    logSpacing = inLogSpacing;

    // Input data
    std::ifstream file(filePath);
    if (!file.is_open())
    {
        std::cout << "File did not open\n";
    }

    std::string linestr; // String for each line to be read into
    std::string data; // String for each data value
    
    // First line will be xVals
    std::getline(file, linestr);

    // Put into string stream, allows us to grab one value at a time
    std::stringstream xStream(linestr);
    
    // Loop through this row
    while (std::getline(xStream, data))
    {
        xVals.push_back(std::stod(data));
    }

    // Next line is probVals. Repeat all steps.
    std::getline(file, linestr);

    std::stringstream probStream(linestr);
    while (std::getline(probStream, data))
    {
        probVals.push_back(std::stod(data));
    }

    file.close();

    N = xVals.size();

    if (logSpacing)
    {
        xMin = log10(xVals[0]);
        xMax = log10(xVals[N-1]);
    }
    else
    {
        xMin = xVals[0];
        xMax = xVals[N-1];
    }

    dx = (xMax - xMin) / (N-1);
}

NumericalDistb::~NumericalDistb()
{

}


double NumericalDistb::evaluateDistribution(double x)
{
    // If a log hist, convert to log space
    if (logSpacing)
    {
        x = log10(x);
    }

    // Throw error if outside bounds
    if ((x < xMin) || (x > xMax))
    {
        std::cout << "Tried to sample distribution outside of bounds\n";
    }

    // Find which index it falls above
    int index = floor((x - xMin)/dx);

    // In rare edge case where x = xMax, clamp downwards
    if (index >= N - 1)
    {
        index = N - 2;
    }

    // Interpolate between the two values
    // (uses math instead of indexing to allow generality for log)
    double x0 = xMin + index*dx;
    double weight = (x - x0) / dx;

    return (probVals[index]*(1.0 - weight) + probVals[index + 1]*(weight));
}