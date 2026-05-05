#pragma once

#include <string>

class ElectronDistb
{
    private:
        double (ElectronDistb::*distbFunc)(double);
        
        // Maxwell boltzmann
        double MBdistb(double beta);
        
        // Maxwell Juttner
        double K1; // Define this outside so it only gets calc'd once
        
        double MJdistb(double beta);
        
    public:
        double T;
        std::string distb; 

        ElectronDistb(double inTheta);

        ~ElectronDistb();

        void updateTemp(double inTheta);

        double operator () (double beta);

        double max();

        void setDistb(int distbNum);
};