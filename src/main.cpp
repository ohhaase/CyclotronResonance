#include "global_vars.hpp"

#include "core/sim_types.hpp"
#include <iostream>
#include <fstream>
#include "nlohmann/json.hpp"

// ===== Main =====
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cout << "Not enough arguments: Missing file path." << std::endl;
        return 1;
    }

    std::string filepath = argv[1];

    std::cout << filepath << std::endl;

    std::ifstream inputFile(filepath);
    nlohmann::json inputParams = nlohmann::json::parse(inputFile);

    Nparticles = inputParams["Nparticles"].get<int>();
    Nbins = inputParams["Nbins"].get<int>();
    simType = inputParams["simType"].get<int>();
    Nthreads = inputParams["Nthreads"].get<int>();

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
        default:
        {
            std::cout << "Error with simtype\n";
        }
    }

    return 0;
}