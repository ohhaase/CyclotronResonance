#pragma once

#include "global_vars.hpp"

#include "Histogram.hpp"
#include "Histogram2D.hpp"
#include "AvgOutput.hpp"

#include <vector>

#include "nlohmann/json.hpp"

enum struct VALTYPE
{
    NRG,
    THETA,
    COUNT,
    POL,
    BETA,
    INIT_NRG,
    INIT_THETA,
    INIT_POL
};

struct HistInfo
{
    Histogram hist;
    std::string name;
    VALTYPE val;
};

struct Hist2DInfo
{
    Histogram2D hist;
    std::string name;
    VALTYPE val_x;
    VALTYPE val_y;
};

struct AvgInfo
{
    AvgOutput avg;
    std::string name;
};

struct HistBounds
{
    double min;
    double max;
};

struct SimData
{
    PhotonState photon;
    PhotonState initPhoton;
    double beta;
};


class OutputHandler
{
    /*
    Class that initializes, stores, and outputs all the data collection in the sim
    - Uses the input file to determine which outputs to create. Can be:
         - Histogram
         - 2D Histogram
         - Avg output
    - Initializes output objects
    - During the simulation, stores values into these objects. Can happen:
         - After each scatter
         - After photon escapes
    - At the end of the simulation, tells all output objects to write their files.
    */

    private:

        std::vector<HistInfo> perScatterHists;
        std::vector<Hist2DInfo> perScatterHists2D;
        
        std::vector<HistInfo> perEscapeHists;
        std::vector<Hist2DInfo> perEscapeHists2D;
        std::vector<AvgInfo> perEscapeAvgs;
        
        void generateHists(nlohmann::json histsInfo);
        void generate2DHists(nlohmann::json hists2DInfo);
        void generateAvgs(nlohmann::json avgsInfo);
        
        HistBounds getBounds(nlohmann::json boundsInfo);
        VALTYPE stringToValtype(std::string val);
        double getValForHist(VALTYPE val, SimData& currentState);
    
    public:
    
        // Constructor (creates objects)
        OutputHandler(nlohmann::json outputParams);


        // Update objects
        void perScatterOutputs(SimData& currentState);

        void perEscapeOutputs(SimData& currentState);


        // Write objects
        void writeOutputs(std::string folder);
};