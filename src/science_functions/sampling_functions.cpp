#include "sampling_functions.hpp"

#include "../global_vars.hpp"
#include "cross_sections.hpp"
#include "../core/helper_functions.hpp"
#include <cmath>


// ===== Sampling functions =====
double sampleDiffCrossSection(double mu_i, double w_i, diffCrossSecFunc diffCrossSec)
{
    // Max probability is at mu_f = 0 or pi
    double maxProb = diffCrossSec(0, mu_i, w_i); 

    // Accept reject
    double x = 1;
    double probability = 0;
    double mu_f;

    while (x > probability)
    {
        mu_f = getRandom(-1, 1);

        probability = diffCrossSec(mu_f, mu_i, w_i);

        x = getRandom(0, maxProb);
    }

    return mu_f;
}

double sampleSigB()
{
    double y = getRandom(cumSigBMin, cumSigBMax);

    // Uses Netwon Raphson method for root finding
    // xn+1 = xn - f(xn)/f'(xn)
    // Func is cumSigmaB - y, derivative is SigmaB
    // TODO: Can we do error in terms of log so that it gives lower values good error too?

    // Initial guess will always be wB
    double xn = wB;

    double xn1 = xn - (cumSigmaB(xn) - y)/SigmaB(xn);

    const double tolerance = 1e-6;
    double error = fabs(xn1 - xn);

    while (error > tolerance)
    {
        xn = xn1;

        xn1 = xn - (cumSigmaB(xn) - y)/SigmaB(xn);

        error = fabs(xn1 - xn);
    }
    
    return xn;
}


bool crossSecAcceptReject(double mu_i, double targetOmega, crossSecFunc crossSec, double max)
{
    double probability = crossSec(mu_i, targetOmega);

    double x = getRandom(0, max);

    if (x < probability)
    {
        return true;
    }
    return false;
}