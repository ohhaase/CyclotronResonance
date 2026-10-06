#include "PhotonDistb.hpp"

#include <iostream>
#include "global_vars.hpp"
#include "core/helper_functions.hpp"

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