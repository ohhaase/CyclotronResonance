#include "global_vars.hpp"

#include "core/sim_types.hpp"
#include <iostream>
#include <fstream>
#include "nlohmann/json.hpp"
#include "helper_objects/OutputHandler.hpp"

// ===== Main =====
int main(int argc, char** argv)
{
    // Confirm correct number of inputs
    if (argc < 2)
    {
        std::cout << "Not enough arguments: Missing file path." << std::endl;
        return 1;
    }

    // Get input json file
    std::string filepath = argv[1];

    std::ifstream inputFile(filepath);
    nlohmann::json inputParams = nlohmann::json::parse(inputFile);

    Nparticles = inputParams["Nparticles"].get<int>();
    Nbins = inputParams["Nbins"].get<int>();
    simType = inputParams["simType"].get<int>();
    Nthreads = inputParams["Nthreads"].get<int>();

    // Initiate sim based on sim type
    switch (simType)
    {
        case 1:
        {
            simType1();
            break;
        }
        case 2:
        {
            simType2();
            break;
        }
        case 3:
        {
            simType3();
            break;
        }
        case 4:
        {
            simType4();
            break;
        }
        case 5:
        {
            std::cout << "Not yet implemented.\n";
            break;
        }
        case 6:
        {
            simType6();
            break;
        }
        case 7:
        {
            // // Temporary test case
            // OutputHandler handler(inputParams["output"]);

            // for (HistInfo info : handler.perScatterHists)
            // {
            //     std::cout << info.hist.minVal << " " << info.hist.maxVal << "\n";
            // }

            break;
        }
        default:
        {
            std::cout << "Error with simtype\n";
        }
    }

    return 0;
}