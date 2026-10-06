#include "PhotonDistb.hpp"

#include <iostream>
#include <string>
#include "nlohmann/json.hpp"
#include "global_vars.hpp"
#include "core/helper_functions.hpp"


PhotonDistb::PhotonDistb(nlohmann::json inputDistbJSON)
{
    // TODO: Maybe some stuff to check that the json is setup right?

    std::string distName = inputDistbJSON["type"].get<std::string>();

    if (distName == "uniform") distbType = PHOTONDISTBFUNC::Uniform;
}


void PhotonDistb::setInputDistribution(PHOTONDISTBFUNC dist)
{
    distbType = dist;
}


double PhotonDistb::sample()
{
    switch (distbType)
    {
        case PHOTONDISTBFUNC::Uniform:
        {
            return uniformSample();
        }

        default:
        {
            std::cout << "Unknown photon distribution type!\n";
        }
    }
}



double PhotonDistb::uniformSample()
{
    return getRandom(lowerOmega, upperOmega);
}