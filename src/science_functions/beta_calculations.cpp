#include "beta_calculations.hpp"

#include "../global_vars.hpp"
#include <cmath>
#include "../helper_objects/ElectronDistb.hpp"
#include "../core/helper_functions.hpp"

// ===== Beta calculations =====
double calcBeta(PhotonState& photon, double w_i, double sign)
{
    double eps_i = photon.omega;
    double Theta_i = photon.theta;

    double Delta = w_i/eps_i;
    double Mu_i = cos(Theta_i);

    double beta = (Mu_i + sign * Delta * sqrt(Delta*Delta + Mu_i*Mu_i - 1)) / (Delta*Delta + Mu_i*Mu_i);

    return beta;
}

double betaForScatter(PhotonState& photon, double w_i)
{
    // Calculate both betas and then decide which one to scatter

    // Calculate plus and minus solutions
    double betaPlus = calcBeta(photon, w_i, 1.0);
    double betaMinus = calcBeta(photon, w_i, -1.0);

    // Calculate probabilities of finding each from distb
    // double probPlus = MBdistb(betaPlus);
    // double probMinus = MBdistb(betaMinus);

    // MJ distb
    double probPlus = electronDistb(betaPlus);
    double probMinus = electronDistb(betaMinus);

    // Choose which one randomly
    double x = getRandom(0, probPlus + probMinus);

    if (x < probPlus)
    {
        return betaPlus;
    }
    else
    {
        return betaMinus;
    }
}

bool radicalConditionMet(PhotonState& photon, double w_i)
{
    double Mu_i = cos(photon.theta);
    double Delta = w_i/photon.omega;

    return (Delta*Delta + Mu_i*Mu_i - 1) >= 0;
}