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
    // Set num of threads
    omp_set_num_threads(Nthreads);

    // Openmp parallel block
    #pragma omp parallel
    {
        // Intialize threadlocal histograms
        Histogram nrgHist{Nbins, lowerOmega, upperOmega};
        Histogram betaHist{Nbins, -1, 1};

        // Loop over each particle we want in the bin
        #pragma omp for
        for (int k = 0; k < Nparticles; k++)
        {
            // Initialize photon
            PhotonState photon{omega, theta, 0, polarization}; 

            // Initialize per-particle histograms
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

            // Once escaped, combine histograms with thread local hist
            nrgHist.combineData(localNRGHist);
            betaHist.combineData(localBetaHist);
        }

        // Combine thread hist with global histograms
        #pragma omp critical
        {
            totalNrgHist.combineData(nrgHist);
            totalBetaHist.combineData(betaHist);
        }
    }
}

AvgPhotonState avgAndBinNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, 
    Histogram& everyNRGHist_in, Histogram& finalNRGHist_in, Histogram& everyThetaHist_in, Histogram& finalThetaHist_in,
    Histogram& everyBetaHist_in, Histogram& finalCountHist_in, Histogram2D& nrgXnrgHist2D_in, Histogram2D& thetaXnrgHist2D_in,
    Histogram2D& nrgXthetaHist2D_in, Histogram2D& thetaXthetaHist2D_in, Histogram2D& finalValsHist2D_in)
{
    // Initialize average results for each
    double avgOmega = 0;
    double avgTheta = 0;
    double avgNum = 0;
    double avgPol = 0;

    omp_set_num_threads(Nthreads);

    // Openmp parallel block
    #pragma omp parallel 
    {
        // Intialize threadlocal histograms
        Histogram everyNRGHist_tl{Nbins, lowerOmega, upperOmega};
        Histogram finalNRGHist_tl{Nbins, finalNRGlow, finalNRGhigh, true};
        Histogram everyThetaHist_tl{Nbins, 0, M_PI};
        Histogram finalThetaHist_tl{Nbins, 0, M_PI};
        Histogram everyBetaHist_tl{Nbins, -1, 1};
        Histogram finalCountHist_tl{Nbins, 1, 100000, true};

        Histogram2D nrgXnrgHist2D_tl{Nbins, Nbins, lowerOmega, upperOmega, finalNRGlow, finalNRGhigh, false, true};
        Histogram2D nrgXthetaHist2D_tl{Nbins, Nbins, lowerOmega, upperOmega, 0, M_PI};
        Histogram2D thetaXnrgHist2D_tl{Nbins, Nbins, 0, M_PI, finalNRGlow, finalNRGhigh, false, true};
        Histogram2D thetaXthetaHist2D_tl{Nbins, Nbins, 0, M_PI, 0, M_PI};
        Histogram2D finalValsHist2D_tl{Nbins, Nbins, finalNRGlow, finalNRGhigh, 0, M_PI, true, false};

        // Loop over each particle we want in the bin
        #pragma for reduction(+: avgOmega, avgTheta, avgNum, avgPol)
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

            // Combine data for per-scatter hists with threadlocal variants
            everyNRGHist_tl.combineData(localNRGHist);
            everyThetaHist_tl.combineData(localThetaHist);
            everyBetaHist_tl.combineData(localBetaHist);

            // Add data for on-escape hists to threadlocal variants
            finalNRGHist_tl.addVal(photon.omega, photon.polarization);
            finalThetaHist_tl.addVal(photon.theta, photon.polarization);
            finalCountHist_tl.addVal(photon.numScatterings, photon.polarization);

            nrgXnrgHist2D_tl.addVal(omega, photon.omega, photon.polarization);
            nrgXthetaHist2D_tl.addVal(omega, photon.theta, photon.polarization);
            thetaXnrgHist2D_tl.addVal(theta, photon.omega, photon.polarization);
            thetaXthetaHist2D_tl.addVal(theta, photon.theta, photon.polarization);

            finalValsHist2D_tl.addVal(photon.omega, photon.theta, photon.polarization);
        }
        
        // Once all particles are done, combine with global hists
        #pragma omp critical // Needs to be in critical block to avoid race conditions
        {
            // Combine data for per-scatter hists
            everyNRGHist_in.combineData(everyNRGHist_tl);
            everyThetaHist_in.combineData(everyThetaHist_tl);
            everyBetaHist_in.combineData(everyBetaHist_tl);

            // Combine data for final value hists
            finalNRGHist_in.combineData(finalNRGHist_tl);
            finalThetaHist_in.combineData(finalThetaHist_tl);
            finalCountHist_in.combineData(finalCountHist_tl);

            nrgXnrgHist2D_in.combineData(nrgXnrgHist2D_tl);
            nrgXthetaHist2D_in.combineData(nrgXthetaHist2D_tl);
            thetaXnrgHist2D_in.combineData(thetaXnrgHist2D_tl);
            thetaXthetaHist2D_in.combineData(thetaXthetaHist2D_tl);

            finalValsHist2D_in.combineData(finalValsHist2D_tl);
        }
    }
    
    // Calculate averages
    avgOmega /= Nparticles;
    avgTheta /= Nparticles;
    avgNum /= Nparticles;
    avgPol /= Nparticles;

    return AvgPhotonState{avgOmega, avgTheta, avgNum, avgPol};
}