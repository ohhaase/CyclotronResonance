#pragma once

enum struct PHOTONDISTBFUNC
{
    Uniform
};


class PhotonDistb
{
    public:

        void setInputDistribution(PHOTONDISTBFUNC dist);

        double sample();

    private:

        PHOTONDISTBFUNC distbType;
        
        
        // Sampling functions

        double uniformSample();

        
};