#pragma once

#include "global_vars.hpp"

#include "Histogram.hpp"
#include "Histogram2D.hpp"
#include "AvgOutput.hpp"

#include <vector>

#include "nlohmann/json.hpp"

struct HistInfo
{
    Histogram hist;
    std::string name;
    std::string val;
};


struct Hist2DInfo
{
    Histogram2D hist;
    std::string name;
    std::string val_x;
    std::string val_y;
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

        std::vector<HistInfo> perScatterHists;
        std::vector<Hist2DInfo> perScatterHists2D;
        
        std::vector<HistInfo> perEscapeHists;
        std::vector<Hist2DInfo> perEscapeHists2D;
        std::vector<AvgInfo> perEscapeAvgs;
        
        void generateHists(nlohmann::json histsInfo);
        void generate2DHists(nlohmann::json hists2DInfo);
        void generateAvgs(nlohmann::json avgsInfo);
        
        HistBounds getBounds(nlohmann::json boundsInfo);
        double getValForHist(std::string val, PhotonState photon, double beta);
    
    public:
    
        // Constructor (creates objects)
        OutputHandler(nlohmann::json outputParams);


        // Update objects
        void perScatterOutputs(PhotonState photon, double beta);

        void perEscapeOutputs(PhotonState initPhoton, PhotonState photon);


        // Write objects
        void writeOutputs();
};