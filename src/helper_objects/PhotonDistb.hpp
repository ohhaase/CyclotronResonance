#pragma once

#include "nlohmann/json.hpp"

enum struct PHOTONDISTBFUNC
{
    Uniform
};


class PhotonDistb
{
    public:
        PhotonDistb(nlohmann::json inputDistbJSON);

        void setInputDistribution(PHOTONDISTBFUNC dist);

        double sample();

    private:

        PHOTONDISTBFUNC distbType;
        
        
        // Sampling functions

        double uniformSample();

        
};