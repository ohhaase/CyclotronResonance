#include "OutputHandler.hpp"

#include "nlohmann/json.hpp"

#include <iostream>
#include <string>

#include "global_vars.hpp"
#include "helper_objects/Histogram.hpp"
#include "helper_objects/Histogram2D.hpp"




OutputHandler::OutputHandler(nlohmann::json outputParams)
{
    // JSON that comes in is expected to be of the form:
    // "hists":
    //      "hist1": <see generateHists() for structure>
    //      "hist2": ...
    //      ...
    // "hists2D":
    //      "hist2D_1": <see generate2DHists() for structure>
    //      "hist2D_2": ...
    //      ...
    // "avgs": <see generateAvgHists() for structure>

    // Possible values of "val" and the types to use for other inputs:
    //      "nrg" (double), "theta" (double), "count" (int), "pol" (double), "beta" (double)

    // Possible values of "trigger":
    //      0: per scatter, 1: per escape

    

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
        generateAvgHists(outputParams["avgs"]);
    }
}


void OutputHandler::generateHists(nlohmann::json histsInfo)
{
    // histsInfo is a json dict where each key is the name of the hist and has the following structure:

    //  "histName": 
    //      "trigger": <int> (see constructor for possible values of "trigger")
    //      "val": <string> (see constructor for possible values of "val")
    //      "Nbins": <int>
    //      "log": <bool>
    //      "bounds": <json> (see getBounds() function)

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

        HistInfo thisHistInfo{thisHist, histIterator.key(), val};

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
    // hists2DInfo is a json dict where each key is the name of the hist2D and has the following structure:

    //  "hist2DName":
    //      "trigger": <int> (see constructor for possible values of "trigger")
    //      "val_x": <string> (see constructor for possible values of "val")
    //      "Nbins_x": <int>
    //      "log_x": <bool>
    //      "bounds_x": <json> (see getBounds() function)
    //      "val_y": <string> (see constructor for possible values of "val")
    //      "Nbins_y": <int>
    //      "log_y": <bool>
    //      "bounds_y": <json> (see getBounds() function)

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

        Hist2DInfo thisHistInfo{thisHist, hist2DIterator.key(), val_x, val_y};

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


void OutputHandler::generateAvgHists(nlohmann::json avgsInfo)
{
    // Average outputs are essentially 2D hists that will simply divide their value by the count at the end
    // x-axis is energy, y-axis is angle
    // They should usually have the same angle and energy ranges as the simulation, though this is not hard coded (needs to be the case in the input file)
    // They will always trigger on photon escape

    // avgsInfo is a json dict with the following structure:

    //  "avgsInfo":
    //      "vals": <list of "val"> (see constructor for possible values of "val")
    //      "Nbins_nrg": <int>
    //      "bounds_nrg": <json> (see getBounds() function)
    //      "Nbins_theta": <int>
    //      "bounds_theta": <json> (see getBounds() function)


    std::vector<std::string> vals = avgsInfo["vals"].get<std::vector<std::string>>();
    int Nbins_nrg = avgsInfo["Nbins_nrg"].get<int>();
    HistBounds bounds_nrg = getBounds(avgsInfo["bounds_nrg"]);
    int Nbins_theta = avgsInfo["Nbins_theta"].get<int>();
    HistBounds bounds_theta = getBounds(avgsInfo["bounds_theta"]);

    // Loop over vals and make histograms
    for (std::string val : vals)
    {
        Histogram2D thisHist{Nbins_nrg, Nbins_theta, bounds_nrg.min, bounds_nrg.max, bounds_theta.min, bounds_theta.max};

        AvgInfo thisAvgInfo{thisHist, val};

        perEscapeAvgs.push_back(thisAvgInfo);
    }
}


HistBounds OutputHandler::getBounds(nlohmann::json boundsInfo)
{
    // Get's the bounds from the input dictionary
    // If a static type, this func references those values from the global state (to be changed)
    // If custom, this func unpacks the values

    // Input JSON has three possible keys:
    // "type": <string> (can be "nrg", "theta", "beta", "custom")
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