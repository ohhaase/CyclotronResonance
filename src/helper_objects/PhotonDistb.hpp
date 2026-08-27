#pragma once

class PhotonDistb
{
    public:

        void setInputDistribution(int dist);

        double sample();

    private:

        int inputType;
        
        // Sampling functions

        double uniformSample();

        
};