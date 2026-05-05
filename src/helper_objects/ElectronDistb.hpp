#pragma once

class ElectronDistb
{
    private:
        double (ElectronDistb::*distbFunc)(double);
        double T;

        // Maxwell boltzmann
        double MBdistb(double beta);

        // Maxwell Juttner
        double K1; // Define this outside so it only gets calc'd once

        double MJdistb(double beta);

    public:
        ElectronDistb(double inTheta);

        ~ElectronDistb();

        void updateTemp(double inTheta);

        double operator () (double beta);

        double max();

        void setDistb(int distbNum);
};