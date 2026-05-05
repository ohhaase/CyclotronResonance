#include "ElectronDistb.hpp"

#define _USE_MATH_DEFINES
#include <cmath>


// Maxwell boltzmann
double ElectronDistb::MBdistb(double beta)
{
    return sqrt(1.0 / (2 * M_PI * T)) * exp(-beta*beta / (2.0 * T));
}


// Maxwell Juttner
double ElectronDistb::MJdistb(double beta)
{
    double gamma = 1.0 / sqrt(1.0 - beta*beta);

    // Uses modified bessel function of 2nd kind
    return (exp(-gamma / T) * gamma*gamma*gamma) / (2.0 * K1);
}


ElectronDistb::ElectronDistb(double inTheta)
{
    // Selects MJ by default
    distbFunc = &ElectronDistb::MJdistb;
    // distbFunc = &ElectronDistb::MBdistb;

    T = inTheta;
    K1 = std::cyl_bessel_k(1.0, 1.0 / T); 
}


ElectronDistb::~ElectronDistb()
{

}


void ElectronDistb::updateTemp(double inT)
{
    T = inT;
    K1 = std::cyl_bessel_k(1.0, 1.0 / T);
}


double ElectronDistb::operator () (double beta)
{
    // Distribtution is callable with (), i.e. electronDistb(beta)
    return (this->*distbFunc)(beta);
}


double ElectronDistb::max()
{
    // Flat out INCORRECT for high Theta. Will need to change.
    return (this->*distbFunc)(0);
}


void ElectronDistb::setDistb(int distbNum)
{
    switch (distbNum)
    {
        case 0: // MB
        {
            distbFunc = &ElectronDistb::MBdistb;
            break;
        }
        case 1: // MJ
        {
            distbFunc = &ElectronDistb::MJdistb;
            break;
        }
    }
}