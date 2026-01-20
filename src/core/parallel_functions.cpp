#include "parallel_functions.hpp"

#include "../global_vars.hpp"
#include "omp.h"
#include "../science_functions/scattering_functions.hpp"
#include "../helper_objects/Histogram.hpp"
#include "../helper_objects/Histogram2D.hpp"
#include "helper_functions.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

// ===== Parallel particle functions =====
AvgPhotonState averageNParticles(double omega, double theta, int polarization, int recoil, int Nparticles)
{
    // Initialize average results for each
    double avgOmega = 0;
    double avgTheta = 0;
    double avgNum = 0;
    double avgPol = 0;

    omp_set_num_threads(Nthreads);

    // Loop over each particle we want in the bin
    #pragma omp parallel for reduction(+: avgOmega, avgTheta, avgNum, avgPol)
    for (int k = 0; k < Nparticles; k++)
    {
        PhotonState photon{omega, theta, 0, polarization}; 

        bool escaped = false; 
        double beta = 0;

        while(!escaped)
        {
            beta = performScatter(photon, recoil);
            escaped = hasEscaped(photon);
        }

        // Once escaped, average results
        avgOmega += photon.omega;
        avgTheta += photon.theta;
        avgNum += photon.numScatterings;
        avgPol += photon.polarization;
    }

    // Calculate averages
    avgOmega /= Nparticles;
    avgTheta /= Nparticles;
    avgNum /= Nparticles;
    avgPol /= Nparticles;

    return AvgPhotonState{avgOmega, avgTheta, avgNum, avgPol};
}

void binNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, Histogram& totalNrgHist, Histogram& totalBetaHist)
{
    // Intialize total histograms
    Histogram nrgHist{Nbins, lowerOmega, upperOmega};
    Histogram betaHist{Nbins, -1, 1};

    omp_set_num_threads(Nthreads);

    // Loop over each particle we want in the bin
    #pragma omp parallel for
    for (int k = 0; k < Nparticles; k++)
    {
        // Initialize photon
        PhotonState photon{omega, theta, 0, polarization}; 

        // Initialize histograms we want
        Histogram localNRGHist{Nbins, lowerOmega, upperOmega};
        Histogram localBetaHist{Nbins, -1, 1};

        // Initialize values that change each loop
        bool escaped = false; 
        double beta = 0;

        // Loop until escaped
        while(!escaped)
        {
            beta = performScatter(photon, recoil);

            // Add values to histograms
            localNRGHist.addVal(photon.omega, photon.polarization);
            localBetaHist.addVal(beta, photon.polarization);

            escaped = hasEscaped(photon);
        }

        // Once escaped, combine histograms
        #pragma omp critical // Needs to be in critical block to avoid race conditions
        {
            nrgHist.combineData(localNRGHist);
            betaHist.combineData(localBetaHist);
        }
    }

    // Combine with global histograms
    totalNrgHist.combineData(nrgHist);
    totalBetaHist.combineData(betaHist);
}

AvgPhotonState avgAndBinNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, 
    Histogram& everyNRGHist, Histogram& finalNRGHist, Histogram& everyThetaHist, Histogram& finalThetaHist,
    Histogram& everyBetaHist, Histogram& finalCountHist, Histogram2D& nrgXnrgHist2D, Histogram2D& thetaXnrgHist2D,
    Histogram2D& nrgXthetaHist2D, Histogram2D& thetaXthetaHist2D, Histogram2D& finalValsHist2D)
{
    // Initialize average results for each
    double avgOmega = 0;
    double avgTheta = 0;
    double avgNum = 0;
    double avgPol = 0;

    omp_set_num_threads(Nthreads);

    // Loop over each particle we want in the bin
    #pragma omp parallel for reduction(+: avgOmega, avgTheta, avgNum, avgPol)
    for (int k = 0; k < Nparticles; k++)
    {
        // Initialize photon state
        int thisPol = polarization;
        if (thisPol == 0)
        {
            double x = getRandom(0, 1);
            if (x < 0.5) thisPol = 1;
        }

        PhotonState photon{omega, theta, 0, thisPol}; 

        // Initialize histograms we want on a per-scatter basis
        Histogram localNRGHist{Nbins, lowerOmega, upperOmega};
        Histogram localThetaHist{Nbins, 0, M_PI};
        Histogram localBetaHist{Nbins, -1, 1};

        bool escaped = false; 
        double beta = 0;

        while(!escaped)
        {
            beta = performScatter(photon, recoil);

            // Add values to histograms
            localNRGHist.addVal(photon.omega, photon.polarization);
            localThetaHist.addVal(photon.theta, photon.polarization);
            localBetaHist.addVal(beta, photon.polarization);

            escaped = hasEscaped(photon);
        }

        // Once escaped, add results to averages
        avgOmega += photon.omega;
        avgTheta += photon.theta;
        avgNum += photon.numScatterings;
        avgPol += photon.polarization;

        #pragma omp critical // Needs to be in critical block to avoid race conditions
        {
            // Combine data for per-scatter hists
            everyNRGHist.combineData(localNRGHist);
            everyThetaHist.combineData(localThetaHist);
            everyBetaHist.combineData(localBetaHist);

            // Add data for on-escape hists
            finalNRGHist.addVal(photon.omega, photon.polarization);
            finalThetaHist.addVal(photon.theta, photon.polarization);
            finalCountHist.addVal(photon.numScatterings, photon.polarization);

            nrgXnrgHist2D.addVal(omega, photon.omega, photon.polarization);
            nrgXthetaHist2D.addVal(omega, photon.theta, photon.polarization);
            thetaXnrgHist2D.addVal(theta, photon.omega, photon.polarization);
            thetaXthetaHist2D.addVal(theta, photon.theta, photon.polarization);

            finalValsHist2D.addVal(photon.omega, photon.theta, photon.polarization);
        }
    }

    // Calculate averages
    avgOmega /= Nparticles;
    avgTheta /= Nparticles;
    avgNum /= Nparticles;
    avgPol /= Nparticles;

    return AvgPhotonState{avgOmega, avgTheta, avgNum, avgPol};
}