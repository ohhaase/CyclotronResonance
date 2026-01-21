#include "scattering_functions.hpp"

#include "../global_vars.hpp"
#include <iostream>
#include <cmath>
#include "cross_sections.hpp"
#include "sampling_functions.hpp"
#include "../core/helper_functions.hpp"
#include "beta_calculations.hpp"

// ===== Scattering Functions =====
bool hasEscaped(PhotonState& photon)
{
    // Prescribed boundaries
    if ((lowerOmega > photon.omega) || (photon.omega > upperOmega)) return true;

    // Make sure things aren't nan
    if (std::isnan(photon.omega) || std::isnan(photon.theta)) return true;

    // if (photon.numScatterings > 100) return true;

    return false;
}


double abberation(double beta, double mu_i)
{
    return (mu_i - beta) / (1.0 - beta * mu_i);
}


PhotonState scatterPhoton(PhotonState& photon, double w_i, double beta, int scatterType, int recoil)
{
    // Determine polarization info
    diffCrossSecFunc diffCrossSec;
    int finalPol;

    switch (scatterType) // Recall: 0 = parallel, 1 = perp
    {
        case 0: // Avg
        {
            diffCrossSec = averageMagDiffCrossSec;
            finalPol = photon.polarization;
            break;
        }
        case 1: // pr2pr
        {
            diffCrossSec = pr2prDiffCrossSec;
            finalPol = 1;
            break;
        }
        case 2: // pr2pl
        {
            diffCrossSec = pr2plDiffCrossSec;
            finalPol = 0;
            break;
        }
        case 3: // pl2pr
        {
            diffCrossSec = pl2prDiffCrossSec;
            finalPol = 1;
            break;
        }
        case 4: //pl2pl
        {
            diffCrossSec = pl2plDiffCrossSec;
            finalPol = 0;
            break;
        }
        default:
        {
            std::cout << "Something went wrong with scatter type!\n";
        }
    }

    // Initial angle in plasma frame
    double Mu_i = cos(photon.theta);

    // Boost angle to ERF
    double mu_i = abberation(beta, Mu_i);

    // Sample final angle
    double mu_f = sampleDiffCrossSection(mu_i, w_i, diffCrossSec);

    // Energy change
    double w_f;
    
    switch (recoil)
    {
        case 0: // No recoil
        {
            w_f = w_i; // Thompson scattering
            break;
        }
        case 1: // Full QED
        {
            if (fabs(mu_f) == 1.0)
            {
                // Special version for singularity
                w_f = 0.5 * ((w_i*w_i * (1.0 - mu_i*mu_i) + 2.0*w_i) / (1.0 + w_i*(1.0 - mu_i)));
            }
            else
            {
                // Otherwise, full QED version
                double dMu = mu_f - mu_i;
                w_f = (1.0 + w_i*(1.0 - mu_i*mu_f) - sqrt(1.0 + 2.0*w_i*mu_f*dMu + w_i*w_i*dMu*dMu)) / (1.0 - mu_f*mu_f);
            }
            break;
        }
    }

    // Boost angle back to plasma frame
    double Mu_f = abberation(-beta, mu_f); // Inverse abberation is just -beta

    // Boost energy back to plasma frame
    double gamma = 1.0 / sqrt(1.0 - beta*beta);
    double eps_f = w_f * gamma * (1.0 + beta * mu_f);

    return PhotonState{eps_f, acos(Mu_f), photon.numScatterings + 1, finalPol};
}

int decideScatterType(PhotonState& photon, double beta, double w_i)
{
    // First, get angle values in ERF
    double mu_i = abberation(beta, cos(photon.theta));

    // Initialize scatter type that will be returned
    // Recall:
    // 0: Avg
    // 1: pr2pr
    // 2: pr2pl
    // 3: pl2pr
    // 4: pl2pl
    int scatterType;

    // First figure out incoming photon polarization
    switch (photon.polarization)
    {
        case -1: // Should only be used when doing polarization averaged
        {
            scatterType = 0;
            break; 
        }
        case 0: // parallel
        {
            // Need to decide between pl2pl and pl2pr

            // Get two probabilities
            double pl2plProb = pl2plTotalCrossSec(mu_i, w_i);
            double pl2prProb = pl2prTotalCrossSec(mu_i, w_i);

            // Generate random number between 0 and sum
            double x = getRandom(0, pl2plProb + pl2prProb);

            // Decide between them
            if (x < pl2plProb)
            {
                scatterType = 4;
            }
            else
            {
                scatterType = 3;
            }

            break;
        }
        case 1: // perp
        {
            // Need to decide between pr2pl and pr2pr

            // Get two probabilities
            double pr2plProb = pr2plTotalCrossSec(mu_i, w_i);
            double pr2prProb = pr2prTotalCrossSec(mu_i, w_i);

            // Generate random number between 0 and sum
            double x = getRandom(0, pr2plProb + pr2prProb);

            // Decide between them
            if (x < pr2plProb)
            {
                scatterType = 2;
            }
            else
            {
                scatterType = 1;
            }

            break;
        }
    }

    return scatterType;
}


double performScatter(PhotonState& photon, int recoil)
{
    // By default, try to scatter at resonance
    double targetOmega = wB; 

    if (radicalConditionMet(photon, targetOmega))
    {
        // If we can scatter at resonance, then do

        // Get beta
        double beta = betaForScatter(photon, targetOmega);

        // Decide final polarization
        int scatterType = decideScatterType(photon, beta, targetOmega);

        // Perform scatter
        photon = scatterPhoton(photon, targetOmega, beta, scatterType, recoil);

        return beta;
    }
    else
    {
        // std::cout<<"tripped"; // For testing purposes
        
        // Otherwise, we need a new target
        // Accept/reject on total cross sections
        bool foundElectron = false;

        // Initialize scatter type
        // Recall:
        // 0: Avg
        // 1: pr2pr
        // 2: pr2pl
        // 3: pl2pr
        // 4: pl2pl
        int scatterType;
        double beta;

        int j = 0;

        while (!foundElectron)
        {
            j += 1;
            // Choose random energy as target, making sure it meets condition
            while (!radicalConditionMet(photon, targetOmega))
            {
                // targetOmega = getRandom(lowerOmega, upperOmega);
                targetOmega = sampleSigB();
            }

            // Calculate beta for that energy
            beta = betaForScatter(photon, targetOmega);

            // Calculate mu_i for that beta
            double mu_i = abberation(beta, cos(photon.theta));

            // Get probabilities from those
            switch (photon.polarization)
            {
                case -1: // Should only be used when doing polarization averaged
                {
                    double avgProbability = crossSecAcceptReject(mu_i, targetOmega, averageMagTotalCrossSec, averageMagTotalCrossSecMax(mu_i));

                    if (avgProbability != 0.0)
                    {
                        foundElectron = true;
                        scatterType = 0;
                    }
                    break;
                }
                case 0: // parallel
                {
                    // Need to decide between pl2pl and pl2pr

                    double pl2plProb = crossSecAcceptReject(mu_i, targetOmega, pl2plTotalCrossSec, pl2plTotalCrossSecMax(mu_i));
                    double pl2prProb = crossSecAcceptReject(mu_i, targetOmega, pl2prTotalCrossSec, pl2prTotalCrossSecMax(mu_i));

                    bool pl2plChosen = (pl2plProb != 0.0);
                    bool pl2prChosen = (pl2prProb != 0.0);

                    // If either individually are found, choose them
                    if (pl2plChosen)
                    {
                        foundElectron = true;
                        scatterType = 4;
                    }

                    if (pl2prChosen)
                    {
                        foundElectron = true;
                        scatterType = 3;
                    }

                    // In rare case that both are accepted, choose randomly between them
                    if (pl2plChosen && pl2prChosen)
                    {
                        foundElectron = true;
                        double x = getRandom(0, pl2plProb + pl2prProb);
                        
                        if (x < pl2plProb)
                        {
                            // pl2pl
                            scatterType = 4;
                        }
                        else
                        {
                            // pl2pr
                            scatterType = 3;
                        }
                    }

                    // Otherwise, nothing is chosen yet
                    break;
                }
                case 1: // perp
                {
                    // Need to decide between pr2pl and pr2pr

                    double pr2plProb = crossSecAcceptReject(mu_i, targetOmega, pr2plTotalCrossSec, pr2plTotalCrossSecMax(mu_i));
                    double pr2prProb = crossSecAcceptReject(mu_i, targetOmega, pr2prTotalCrossSec, pr2prTotalCrossSecMax(mu_i));

                    bool pr2plChosen = (pr2plProb != 0.0);
                    bool pr2prChosen = (pr2prProb != 0.0);

                    // If either individually are found, choose them
                    if (pr2plChosen)
                    {
                        foundElectron = true;
                        scatterType = 2;
                    }

                    if (pr2prChosen)
                    {
                        foundElectron = true;
                        scatterType = 1;
                    }

                    // In rare case that both are accepted, choose randomly between them
                    if (pr2plChosen && pr2prChosen)
                    {
                        foundElectron = true;
                        double x = getRandom(0, pr2plProb + pr2prProb);
                        
                        if (x < pr2plProb)
                        {
                            // pr2pl
                            scatterType = 2;
                        }
                        else
                        {
                            // pr2pr
                            scatterType = 1;
                        }
                    }

                    // Otherwise, nothing is chosen yet
                    break;
                }
            }
        }
        
        // Once we have a target that satisfies radical and was accepted, scatter using those values
        photon = scatterPhoton(photon, targetOmega, beta, scatterType, recoil);

        // std::cout << j << "\n";

        return beta;
    }

}