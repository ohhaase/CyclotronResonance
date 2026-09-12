#include "OutputHandler.hpp"

#include "nlohmann/json.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>

#include "global_vars.hpp"
#include "helper_objects/Histogram.hpp"
#include "helper_objects/Histogram2D.hpp"
#include "helper_objects/AvgOutput.hpp"




OutputHandler::OutputHandler(nlohmann::json outputParams)
{
    /*
    JSON that comes in is expected to be of the form:
    "hists":
         "hist1": <see generateHists() for structure>
         "hist2": ...
         ...
    "hists2D":
         "hist2D_1": <see generate2DHists() for structure>
         "hist2D_2": ...
         ...
    "avgs": 
         "avg1": <see generateAvgHists() for structure>
         ...

    Possible values of "trigger":
         0: per scatter, 1: per escape
    */
    

    // Create output objects
    if (outputParams.contains("hists"))
    {
        generateHists(outputParams["hists"]);
    }

    if (outputParams.contains("hists2D"))
    {
        generate2DHists(outputParams["hists2D"]);
    }

    if (outputParams.contains("avgs"))
    {
        generateAvgs(outputParams["avgs"]);
    }
}


void OutputHandler::generateHists(nlohmann::json histsInfo)
{
    /*
    histsInfo is a json dict where each key is the name of the hist and has the following structure:

     "histName": 
         "trigger": <int> (see constructor for possible values of "trigger")
         "val": <string> (see stringToValtype() for possible values)
         "Nbins": <int>
         "log": <bool>
         "bounds": <json> (see getBounds() function)
    */

    // Loop over all hist objects
    for (nlohmann::json::iterator histIterator = histsInfo.begin(); histIterator != histsInfo.end(); ++histIterator)
    {
        nlohmann::json thisHistJSON = histIterator.value();

        // Get params
        int trigger = thisHistJSON["trigger"].get<int>();
        std::string val = thisHistJSON["val"].get<std::string>();
        int Nbins = thisHistJSON["Nbins"].get<int>();
        bool isLog = thisHistJSON["log"].get<bool>();
        HistBounds bounds = getBounds(thisHistJSON["bounds"]);

        Histogram thisHist{Nbins, bounds.min, bounds.max, isLog};

        HistInfo thisHistInfo{thisHist, histIterator.key(), stringToValtype(val)};

        switch (trigger)
        {
            case 0:
            {
                // Per scatter
                perScatterHists.push_back(thisHistInfo);
                break;
            }
            case 1:
            {
                // Per escape
                perEscapeHists.push_back(thisHistInfo);
                break;
            }
            default:
            {
                std::cout << "Unknown histogram trigger!\n";
            }
        }
    }
}


void OutputHandler::generate2DHists(nlohmann::json hists2DInfo)
{
    /*
    hists2DInfo is a json dict where each key is the name of the hist2D and has the following structure:

     "hist2DName":
         "trigger": <int> (see constructor for possible values of "trigger")
         "val_x": <string> (see stringToValtype() for possible values)
         "Nbins_x": <int>
         "log_x": <bool>
         "bounds_x": <json> (see getBounds() function)
         "val_y": <string> (see stringToValtype() for possible values)
         "Nbins_y": <int>
         "log_y": <bool>
         "bounds_y": <json> (see getBounds() function)
    */

    // Loop over all hist objects
    for (nlohmann::json::iterator hist2DIterator = hists2DInfo.begin(); hist2DIterator != hists2DInfo.end(); ++hist2DIterator)
    {
        nlohmann::json thisHistJSON = hist2DIterator.value();

        int trigger = thisHistJSON["trigger"].get<int>();

        // Get x params
        int Nbins_x = thisHistJSON["Nbins_x"].get<int>();
        std::string val_x = thisHistJSON["val_x"].get<std::string>();
        bool isLog_x = thisHistJSON["log_x"].get<bool>();
        HistBounds bounds_x = getBounds(thisHistJSON["bounds_x"]);


        // Get y params
        int Nbins_y = thisHistJSON["Nbins_y"].get<int>();
        std::string val_y = thisHistJSON["val_y"].get<std::string>();
        bool isLog_y = thisHistJSON["log_y"].get<bool>();
        HistBounds bounds_y = getBounds(thisHistJSON["bounds_y"]);


        Histogram2D thisHist{Nbins_x, Nbins_y, bounds_x.min, bounds_x.max, bounds_y.min, bounds_y.max, isLog_x, isLog_y};

        Hist2DInfo thisHistInfo{thisHist, hist2DIterator.key(), stringToValtype(val_x), stringToValtype(val_y)};

        switch (trigger)
        {
            case 0:
            {
                // Per scatter
                perScatterHists2D.push_back(thisHistInfo);
                break;
            }
            case 1:
            {
                // Per escape
                perEscapeHists2D.push_back(thisHistInfo);
                break;
            }
            default:
            {
                std::cout << "Unknown histogram trigger!\n";
            }
        }
    }
}


void OutputHandler::generateAvgs(nlohmann::json avgsInfo)
{
    /*
    Average outputs are essentially 2D hists that will simply divide their value by the count at the end
    x-axis is energy, y-axis is angle
    They should usually have the same angle and energy ranges as the simulation, though this is not hard coded (they need to be specified in the input file)
    They will always trigger on photon escape

    avgsInfo is a json dict where each key is the name of the avg output and has the following structure:

     "avg1name":
         "Nbins_nrg": <int>
         "bounds_nrg": <json> (see getBounds() function)
         "Nbins_theta": <int>
         "bounds_theta": <json> (see getBounds() function)
    */

    for (nlohmann::json::iterator avgIterator = avgsInfo.begin(); avgIterator != avgsInfo.end(); ++avgIterator)
    {
        nlohmann::json thisAvgJSON = avgIterator.value();

        int Nbins_nrg = thisAvgJSON["Nbins_nrg"].get<int>();
        HistBounds bounds_nrg = getBounds(thisAvgJSON["bounds_nrg"]);
        int Nbins_theta = thisAvgJSON["Nbins_theta"].get<int>();
        HistBounds bounds_theta = getBounds(thisAvgJSON["bounds_theta"]);

        AvgOutput thisAvg{Nbins_nrg, Nbins_theta, bounds_nrg.min, bounds_nrg.max, bounds_theta.min, bounds_theta.max};

        AvgInfo thisAvgInfo{thisAvg, avgIterator.key()};

        perEscapeAvgs.push_back(thisAvgInfo);
    }
}


HistBounds OutputHandler::getBounds(nlohmann::json boundsInfo)
{
    // Get's the bounds from the input dictionary
    // If a static type, this func references those values from the global state (to be changed)
    // If custom, this func unpacks the values

    // Input JSON has three possible keys:
    // "type": <string> (can be "nrg", "esc_nrg", "theta", "beta", "custom")
    // "min": <double> (min value, only needed if "type"="custom")
    // "max": <double> (max value, only needed if "type"="custom")

    HistBounds bounds{0, 0};

    std::string type = boundsInfo["type"].get<std::string>();

    if (type == "nrg")
    {
        bounds.min = lowerOmega;
        bounds.max = upperOmega;

        return bounds;
    }

    if (type == "esc_nrg")
    {
        bounds.min = finalNRGlow;
        bounds.max = finalNRGhigh;

        return bounds;
    }

    if (type == "theta")
    {
        bounds.min = 0;
        bounds.max = M_PI;

        return bounds;
    }

    if (type == "beta")
    {
        bounds.min = -1;
        bounds.max = 1;

        return bounds;
    }

    if (type == "custom")
    {
        if (boundsInfo.contains("min") & boundsInfo.contains("max"))
        {
            bounds.min = boundsInfo["min"].get<double>();
            bounds.max = boundsInfo["max"].get<double>();
            
            return bounds;
        }
    }

    // If we reach this point, our bounds are improperly defined
    std::cout << "Hist bins improperly defined!\n";

    return bounds; // Better error handling eventually?
}


VALTYPE OutputHandler::stringToValtype(std::string val)
{
    /*
    We include this separate function so that writing values to histograms can be optimized
    Possible values of "val":
         "nrg", "theta", "count", "pol", "beta", "init_nrg", "init_theta", "init_pol"
    */

    if (val == "nrg") return VALTYPE::NRG;

    if (val == "theta") return VALTYPE::THETA;

    if (val == "count") return VALTYPE::COUNT;

    if (val == "pol") return VALTYPE::POL;

    if (val == "beta") return VALTYPE::BETA;

    if (val == "init_nrg") return VALTYPE::INIT_NRG;

    if (val == "init_theta") return VALTYPE::INIT_THETA;

    if (val == "init_pol") return VALTYPE::INIT_POL;

    // If we get to this point, something is awry
    std::cout << "Prescribed val for hist is invalid.\n";
    return VALTYPE::NRG;
} 



double OutputHandler::getValForHist(VALTYPE val, SimData& data)
{        
    switch (val)
    {
        case VALTYPE::NRG:
        {
            return data.photon.omega;
        }
        case VALTYPE::THETA:
        {
            return data.photon.theta;
        }
        case VALTYPE::COUNT:
        {
            return static_cast<double>(data.photon.numScatterings);
        }
        case VALTYPE::POL:
        {
            return static_cast<double>(data.photon.polarization);
        }
        case VALTYPE::BETA:
        {
            return data.beta;
        }
        case VALTYPE::INIT_NRG:
        {
            return data.initPhoton.omega;
        }
        case VALTYPE::INIT_THETA:
        {
            return data.initPhoton.theta;
        }
        case VALTYPE::INIT_POL:
        {
            return data.initPhoton.polarization;
        }
        default:
        {
            std::cout << "Tried to write unknown value to hist\n";
            return -1;
        }
    }
}


void OutputHandler::perScatterOutputs(SimData& data)
{
    // Add all histogram data
    for (HistInfo thisHist : perScatterHists)
    {
        thisHist.hist.addVal(getValForHist(thisHist.val, data), data.photon.polarization);
    }

    // Add all hist2D data
    for (Hist2DInfo thisHist : perScatterHists2D)
    {
        thisHist.hist.addVal(getValForHist(thisHist.val_x, data), getValForHist(thisHist.val_y, data), data.photon.polarization);
    }
}


void OutputHandler::perEscapeOutputs(SimData& data)
{
    // Add all histogram data
    for (HistInfo thisHist : perEscapeHists)
    {
        thisHist.hist.addVal(getValForHist(thisHist.val, data), data.photon.polarization);
    }

    // Add all hist2D data
    for (Hist2DInfo thisHist : perEscapeHists2D)
    {
        thisHist.hist.addVal(getValForHist(thisHist.val_x, data), getValForHist(thisHist.val_y, data), data.photon.polarization);
    }

    // Add all average data
    for (AvgInfo thisAvg : perEscapeAvgs)
    {
        thisAvg.avg.addVal(data.initPhoton, data.photon);
    }
}


void OutputHandler::writeOutputs(std::string folder)
{
    // Write all per scatter hists
    for (HistInfo thisHist : perScatterHists)
    {
        thisHist.hist.exportToFile(thisHist.name, folder);
    }

    // Write all per scatter 2D hists
    for (Hist2DInfo thisHist2D : perScatterHists2D)
    {
        thisHist2D.hist.exportToFile(thisHist2D.name, folder);
    }

    // Write all per escape hists
    for (HistInfo thisHist : perScatterHists)
    {
        thisHist.hist.exportToFile(thisHist.name, folder);
    }

    // Write all per escape 2D hists
    for (Hist2DInfo thisHist2D : perEscapeHists2D)
    {
        thisHist2D.hist.exportToFile(thisHist2D.name, folder);
    }

    // Write all per escape avgs
    for (AvgInfo thisAvg : perEscapeAvgs)
    {
        thisAvg.avg.exportToFile(thisAvg.name, folder);
    }
}