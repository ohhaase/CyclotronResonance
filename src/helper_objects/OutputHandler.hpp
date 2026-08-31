#pragma once

#include "global_vars.hpp"

#include "Histogram.hpp"
#include "Histogram2D.hpp"

#include <vector>

#include "nlohmann/json.hpp"

class OutputHandler
{
    // Class that initializes, stores, and outputs all the data collection in the sim
    // - Uses the input file to determine which outputs to create. Can be:
    //      - Histogram
    //      - Avg output
    // - Initializes output objects
    // - During the simulation, stores values into these objects. Can happen:
    //      - After each scatter
    //      - After photon escapes
    // - At the end of the simulation, tells all output objects to write their files.

    private:

        std::vector<Histogram> perScatterHists;
        std::vector<Histogram2D> perScaterHists2D;

        std::vector<Histogram> perEscapeHists;
        std::vector<Histogram2D> perEscapeHists2D;
        // std::vector<HistogramAvg> perEscapeAvgs;

    public:

        // Constructor (creates objects)
        OutputHandler(nlohmann::json outputParams);


        // Update objects
        void perScatterOutputs();

        void perEscapeOutputs();


        // Write objects
        void writeOutputs();
};