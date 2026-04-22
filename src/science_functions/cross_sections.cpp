#include "cross_sections.hpp"

#include "../global_vars.hpp"
#include <cmath>

extern const double SigmaBMax = 2.0 * wB * wB / (Gamma * Gamma);
double SigmaB(double w)
{
    double wPlus = w + wB;
    double wMinus = w - wB;

    return 0.5 * w*w * (1.0 / (wMinus*wMinus + 0.25*Gamma*Gamma) + 1.0 / (wPlus * wPlus));
}


double cumSigmaB(double w)
{
    double wPlus = wB + w;
    double wMinus = wB - w;

    double A = (Gamma * Gamma - 4 * wB * wB) / (4 * Gamma);
    double B = A * (atan(2.0 * wMinus / Gamma) - atan(2.0 * wPlus / Gamma));
    double C = 0.5 * wB * log((Gamma*Gamma + 4*wMinus*wMinus) / (Gamma*Gamma + 4*wPlus*wPlus));

    return w + B + C;
}
double cumSigBMin = cumSigmaB(lowerOmega);
double cumSigBMax = cumSigmaB(upperOmega);


// -- Polarization averaged (0) --
double averageMagTotalCrossSec(double mu_i, double w_i)
{
    // Given in units of sigma_T
    double sigB = SigmaB(w_i);

    return 0.25 * (1.0 - mu_i*mu_i + sigB*(1 + mu_i*mu_i));
}

double averageMagTotalCrossSecMax(double mu_i)
{
    // Given in units of sigma_T
    double sigB = SigmaBMax;

    return 0.25 * (1.0 - mu_i*mu_i + sigB*(1 + mu_i*mu_i));
}

double averageMagDiffCrossSec(double mu_f, double mu_i, double w_i)
{
    double sigB = SigmaB(w_i);

    return 0.25 * (1.0 + 1.5*mu_f*mu_f + (6.0*(1-mu_f*mu_f)*(1-mu_i*mu_i) 
    + 3.0 * mu_i*mu_i * mu_f*mu_f * sigB)/(4.0 + sigB + (2.0*mu_i*mu_i - 1.0)*(sigB - 4.0)));
}


// -- Perpendicular to Perpendicular (1) --
double pr2prTotalCrossSec(double mu_i, double w_i)
{
    // Given in units of sigma_T
    double sigB = SigmaB(w_i);

    return 0.75 * sigB;
}

double pr2prTotalCrossSecMax(double mu_i)
{
    // Given in units of sigma_T
    double sigB = SigmaBMax;

    return 0.75 * sigB;
}

double pr2prDiffCrossSec(double mu_f, double mu_i, double w_i)
{
    return 0.5;
}


// -- Perpendicular to Parallel (2) --
double pr2plTotalCrossSec(double mu_i, double w_i)
{
    // Given in units of sigma_T
    double sigB = SigmaB(w_i);

    return 0.25 * sigB;
}

double pr2plTotalCrossSecMax(double mu_i)
{
    // Given in units of sigma_T
    double sigB = SigmaBMax;

    return 0.25 * sigB;
}

double pr2plDiffCrossSec(double mu_f, double mu_i, double w_i)
{
    return 1.5 * mu_f * mu_f;
}


// -- Parallel to Perpendicular (3) --
double pl2prTotalCrossSec(double mu_i, double w_i)
{
    // Given in units of sigma_T
    double sigB = SigmaB(w_i);

    return 0.75 * mu_i * mu_i * sigB;
}

double pl2prTotalCrossSecMax(double mu_i)
{
    // Given in units of sigma_T
    double sigB = SigmaBMax;

    return 0.75 * mu_i * mu_i * sigB;
}

double pl2prDiffCrossSec(double mu_f, double mu_i, double w_i)
{
    return 0.5;
}


// -- Parallel to Parallel (4) --
double pl2plTotalCrossSec(double mu_i, double w_i)
{
    // Given in units of sigma_T
    double sigB = SigmaB(w_i);

    return (1.0 - mu_i*mu_i) + 0.25*mu_i*mu_i*sigB;
}

double pl2plTotalCrossSecMax(double mu_i)
{
    // Given in units of sigma_T
    double sigB = SigmaBMax;

    return (1.0 - mu_i*mu_i) + 0.25*mu_i*mu_i*sigB;
}

double pl2plDiffCrossSec(double mu_f, double mu_i, double w_i)
{
    double sigB = SigmaB(w_i);

    double numerator = (6.0 * (1.0 - mu_f*mu_f) * (1.0 - mu_i*mu_i) + 3.0*mu_i*mu_i*mu_f*mu_f*sigB);
    
    double denominator = (4.0 + (2.0 * mu_i*mu_i - 1.0)*(sigB - 4.0) + sigB);

    return numerator / denominator;
}