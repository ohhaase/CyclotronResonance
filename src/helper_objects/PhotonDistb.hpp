#pragma once

enum struct PHOTONDISTBFUNC
{
    Uniform
};


class PhotonDistb
{
    public:

        void setInputDistribution(int dist);

        double sample();

    private:

        PHOTONDISTBFUNC distbType;
        
        
        // Sampling functions

        double uniformSample();

        
};