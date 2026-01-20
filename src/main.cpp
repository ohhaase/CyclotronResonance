#include "global_vars.hpp"

#include "core/sim_types.hpp"
#include <iostream>
#include <fstream>

// ===== Main =====
int main(int argc, char** argv)
{
    if (argc < 4)
    {
        std::cout << "Not enough arguments: [Nparticles, Nbins, simType, Nthreads]" << std::endl;
        return 1;
    }
    else
    {
        Nparticles = std::atoi(argv[1]);
        Nbins = std::atoi(argv[2]);
        simType = std::atoi(argv[3]);
        Nthreads = std::atoi(argv[4]);
    }
    
    std::ofstream simParams("simParams.csv");

    // Write simulation data to file
    simParams << Nparticles << "," << Nbins << "," << simType;

    simParams.close();

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
        default:
        {
            std::cout << "Error with simtype\n";
        }
    }

    return 0;
}