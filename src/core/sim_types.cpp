#include "sim_types.hpp"

#include "../global_vars.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include "../helper_objects/Timer.hpp"
#include "parallel_functions.hpp"
#include "helper_functions.hpp"
#include "../helper_objects/ElectronDistb.hpp"

#define _USE_MATH_DEFINES
#include <cmath>


// ===== Sim functions =====
void simType1()
{
    std::ofstream outputFile("outdata.csv");

    Timer totalTime;

    // Loop over each energy bin
    for (int i = 0; i < Nbins; i++)
    {
        // Linear spacing between the two limits
        double omega = lowerOmega + (double)i/Nbins * (upperOmega - lowerOmega);
        
        // Loop over each angle bin
        for (int j = 0; j < Nbins; j++)
        {
            // Linear spacing from 0 to pi (make it the center of cells)
            double theta = (double)j/Nbins * M_PI;

            AvgPhotonState avgPhoton = averageNParticles(omega, theta, 0, Nparticles, 1);

            // Write these averages to output file with the current hist coordinates
            outputFile << omega << "," << theta << "," << avgPhoton.omega << "," 
                << avgPhoton.theta << "," << avgPhoton.numScatterings << "," << avgPhoton.polarization << "\n";

            loadingBar(Nbins, i, j);
        }
    }
    
    outputFile.close();

    std::cout << "Simulation completed in: " << totalTime.elapsedSec() << " seconds\n";

    system("py dataprocessing.py");
}

void simType2()
{
    // Do we still want to do this looping over the bins?
    // Initialize histograms
    Histogram globalNRGHist{Nbins, lowerOmega, upperOmega};
    Histogram globalBetaHist{Nbins, -1, 1};

    // Loop over each energy bin
    for (int i = 0; i < Nbins; i++)
    {
        // Linear spacing between the two limits
        double omega = lowerOmega + (double)i/Nbins * (upperOmega - lowerOmega);
        
        // Loop over each angle bin
        for (int j = 0; j < Nbins; j++)
        {
            // Linear spacing from 0 to pi
            double theta = (double)j/Nbins * M_PI;

            binNParticles(omega, theta, -1, Nparticles, 1, globalNRGHist, globalBetaHist);

            loadingBar(Nbins, i, j);
        }
    }
    
    // Write histograms to file
    globalNRGHist.exportToFile("energy");
    globalBetaHist.exportToFile("beta");

    system("py dataprocessing.py");
}

void simType3()
{
    // A third sim type to pretty much just time things
    // Will just output timing data:
    // Total sim time (ms)
    // Hist of bin times

    Histogram binTimings{50, 0, 10};

    Timer simTime;

    double avgBinTime = 0;

    // Loop over bins
    for (int i = 0; i < Nbins; i++)
    {
        // Linear spacing between the two limits (at center of bins)
        double omega = lowerOmega + (double)(i+0.5)/Nbins * (upperOmega - lowerOmega);

        for (int j = 0; j < Nbins; j++)
        {
            Timer binTime;

            // Linear spacing from 0 to pi (at center of bins)
            double theta = (double)(j + 0.5)/Nbins * M_PI;

            int pol = 0;
            if (getRandom(0, 1) < 0.5) pol = 1;

            // Bin N particles
            AvgPhotonState avgPhoton = averageNParticles(omega, theta, pol, 1, Nparticles);

            avgBinTime += binTime.elapsedMilli();
            binTimings.addVal(binTime.elapsedMilli(), 0);
        }
    }
    
    avgBinTime /= (Nbins*Nbins);

    // Write these couple numbers to console
    // Comment out top line when you want to port into a csv file
    // std::cout << "Sim finished with [Nparticles, Nbins, Nthreads, Time elapsed (s), Average bin time (ms)] =\n";
    std::cout << Nparticles << "," << Nbins << "," << Nthreads << "," << simTime.elapsedSec() << "," << avgBinTime << "\n";

    binTimings.exportToFile("binTimings");
}


void simType4()
{
    // Outputs both average and histogram data

    // Performs 4 simulations:
    // Theta = 0.05, recoil on
    // Theta = 0.05, recoil off
    // Theta = 0.025, recoil on
    // Theta = 0.025, recoil off

    double Thetas[2] = {0.05, 0.025};
    // double Thetas[2] = {0.01, 0.005};
    int recoils[2] = {0, 1}; 

    // Loop over electron temeperatures
    for (int i = 0; i < 2; i++)
    {
        electronDistb.updateTheta(Thetas[i]);

        // Loop over recoil
        for (int j = 0; j < 2; j++)
        {
            // Unique folder
            std::string folderName = std::to_string(i) + std::to_string(j);

            std::filesystem::create_directory(folderName);

            // Initialize output file for average data
            std::ofstream outputFile(folderName + "/" + "avgdata.csv");

            // Initialize histograms
            Histogram everyNRGHist{Nbins, lowerOmega, upperOmega};
            Histogram finalNRGHist{Nbins, finalNRGlow, finalNRGhigh, true};
            Histogram everyThetaHist{Nbins, 0, M_PI};
            Histogram finalThetaHist{Nbins, 0, M_PI};
            Histogram everyBetaHist{Nbins, -1, 1};
            Histogram finalCountHist{Nbins, 1, 10000, true};

            // Initialize 2D histograms
            Histogram2D nrgXnrgHist2D{Nbins, Nbins, lowerOmega, upperOmega, finalNRGlow, finalNRGhigh, false, true};
            Histogram2D nrgXthetaHist2D{Nbins, Nbins, lowerOmega, upperOmega, 0, M_PI};
            Histogram2D thetaXnrgHist2D{Nbins, Nbins, 0, M_PI, finalNRGlow, finalNRGhigh, false, true};
            Histogram2D thetaXthetaHist2D{Nbins, Nbins, 0, M_PI, 0, M_PI};
            Histogram2D finalValsHist2D{Nbins, Nbins, finalNRGlow, finalNRGhigh, 0, M_PI, true, false};

            // Loop over each energy bin
            for (int ii = 0; ii < Nbins; ii++)
            {
                // Linear spacing between the two limits (at center of bins)
                double omega = lowerOmega + (double)(ii+0.5)/Nbins * (upperOmega - lowerOmega);
                
                // Loop over each angle bin
                for (int jj = 0; jj < Nbins; jj++)
                {
                    // Linear spacing from 0 to pi (at center of bins)
                    double theta = (double)(jj + 0.5)/Nbins * M_PI;

                    // Initialize timer for bin
                    Timer binTimer;

                    // Bin N particles
                    AvgPhotonState avgPhoton = avgAndBinNParticles(omega, theta, 0, Nparticles, recoils[j], 
                        everyNRGHist, finalNRGHist, everyThetaHist, finalThetaHist, everyBetaHist, finalCountHist,
                        nrgXnrgHist2D, thetaXnrgHist2D, nrgXthetaHist2D, thetaXthetaHist2D, finalValsHist2D);

                    // Write these averages to output file with the current hist coordinates
                    outputFile << omega << "," << theta << "," << avgPhoton.omega << "," 
                        << avgPhoton.theta << "," << avgPhoton.numScatterings << "," << avgPhoton.polarization 
                        << "," << binTimer.elapsedMicro() << "\n";

                    loadingBar(Nbins, ii, jj);
                }
            }

            // Close avg file
            outputFile.close();

            // Export histograms
            everyNRGHist.exportToFile("nrg", folderName);
            everyThetaHist.exportToFile("theta", folderName);
            everyBetaHist.exportToFile("beta", folderName);

            finalNRGHist.exportToFile("esc_nrg", folderName);
            finalThetaHist.exportToFile("esc_theta", folderName);
            finalCountHist.exportToFile("num", folderName);

            nrgXnrgHist2D.exportToFile("nrgXnrg", folderName);
            nrgXthetaHist2D.exportToFile("nrgXtheta", folderName);
            thetaXnrgHist2D.exportToFile("thetaXnrg", folderName);
            thetaXthetaHist2D.exportToFile("thetaXtheta", folderName);

            finalValsHist2D.exportToFile("finalVals", folderName);
        }
    }
}


void simType5()
{
    // NOT YET IMPLEMENTED
    // Uses an input distribution, outputs histograms

    // Performs 4 simulations:
    // Theta = 0.05, recoil on
    // Theta = 0.05, recoil off
    // Theta = 0.025, recoil on
    // Theta = 0.025, recoil off

    double Thetas[2] = {0.05, 0.025};
    int recoils[2] = {0, 1}; 

    // Loop over electron temeperatures
    for (int i = 0; i < 2; i++)
    {
        electronDistb.updateTheta(Thetas[i]);

        // Loop over recoil
        for (int j = 0; j < 2; j++)
        {
            // Unique folder
            std::string folderName = std::to_string(i) + std::to_string(j);

            std::filesystem::create_directory(folderName);

            // Initialize histograms
            Histogram everyNRGHist{Nbins, lowerOmega, upperOmega};
            Histogram finalNRGHist{Nbins, finalNRGlow, finalNRGhigh, true};
            Histogram everyThetaHist{Nbins, 0, M_PI};
            Histogram finalThetaHist{Nbins, 0, M_PI};
            Histogram everyBetaHist{Nbins, -1, 1};
            Histogram finalCountHist{Nbins, 0, maxNum};

            // Initialize 2D histogram vectors
            Histogram2D nrgXnrgHist2D{Nbins, Nbins, lowerOmega, upperOmega, finalNRGlow, finalNRGhigh, false, true};
            Histogram2D nrgXthetaHist2D{Nbins, Nbins, lowerOmega, upperOmega, 0, M_PI};
            Histogram2D thetaXnrgHist2D{Nbins, Nbins, 0, M_PI, finalNRGlow, finalNRGhigh, false, true};
            Histogram2D thetaXthetaHist2D{Nbins, Nbins, 0, M_PI, 0, M_PI};
            Histogram2D finalValsHist2D{Nbins, Nbins, finalNRGlow, finalNRGhigh, 0, M_PI, true, false};

            // Get initial energy and angle
            // double omega = sampleOmega();
            // double theta = sampleTheta();

            double omega;
            double theta;

            // Bin N particles - should I switch to just binning? does avging make that much diff?
            AvgPhotonState avgPhoton = avgAndBinNParticles(omega, theta, 0, Nparticles, recoils[j], 
                everyNRGHist, finalNRGHist, everyThetaHist, finalThetaHist, everyBetaHist, finalCountHist,
                nrgXnrgHist2D, thetaXnrgHist2D, nrgXthetaHist2D, thetaXthetaHist2D, finalValsHist2D);

            // Export histograms
            everyNRGHist.exportToFile("nrg", folderName);
            everyThetaHist.exportToFile("theta", folderName);
            everyBetaHist.exportToFile("beta", folderName);

            finalNRGHist.exportToFile("esc_nrg", folderName);
            finalThetaHist.exportToFile("esc_theta", folderName);
            finalCountHist.exportToFile("num", folderName);

            nrgXnrgHist2D.exportToFile("nrgXnrg", folderName);
            nrgXthetaHist2D.exportToFile("nrgXtheta", folderName);
            thetaXnrgHist2D.exportToFile("thetaXnrg", folderName);
            thetaXthetaHist2D.exportToFile("thetaXtheta", folderName);

            finalValsHist2D.exportToFile("finalVals", folderName);
        }
    }
}


void simType6()
{
    // Outputs both average and histogram data
    // Performs 4 different temp simulations (low temp, uses MB)
    // No recoil
    // Intended for pairing with the scatter num cutoffs

    double Thetas[4] = {0.01, 0.007, 0.005, 0.001};

    electronDistb.setDistb(0); // MB distb

    // Loop over electron temeperatures
    for (int i = 0; i < 4; i++)
    {
        electronDistb.updateTheta(Thetas[i]);

        // Unique folder
        std::string folderName = std::to_string(i);

        std::filesystem::create_directory(folderName);

        // Initialize output file for average data
        std::ofstream outputFile(folderName + "/" + "avgdata.csv");

        // Initialize histograms
        Histogram everyNRGHist{Nbins, lowerOmega, upperOmega};
        Histogram finalNRGHist{Nbins, finalNRGlow, finalNRGhigh, true};
        Histogram everyThetaHist{Nbins, 0, M_PI};
        Histogram finalThetaHist{Nbins, 0, M_PI};
        Histogram everyBetaHist{Nbins, -1, 1};
        Histogram finalCountHist{Nbins, 1, 100000, true};

        // Initialize 2D histograms
        Histogram2D nrgXnrgHist2D{Nbins, Nbins, lowerOmega, upperOmega, finalNRGlow, finalNRGhigh, false, true};
        Histogram2D nrgXthetaHist2D{Nbins, Nbins, lowerOmega, upperOmega, 0, M_PI};
        Histogram2D thetaXnrgHist2D{Nbins, Nbins, 0, M_PI, finalNRGlow, finalNRGhigh, false, true};
        Histogram2D thetaXthetaHist2D{Nbins, Nbins, 0, M_PI, 0, M_PI};
        Histogram2D finalValsHist2D{Nbins, Nbins, finalNRGlow, finalNRGhigh, 0, M_PI, true, false};

        // Loop over each energy bin
        for (int ii = 0; ii < Nbins; ii++)
        {
            // Linear spacing between the two limits (at center of bins)
            double omega = lowerOmega + (double)(ii+0.5)/Nbins * (upperOmega - lowerOmega);
            
            // Loop over each angle bin
            for (int jj = 0; jj < Nbins; jj++)
            {
                // Linear spacing from 0 to pi (at center of bins)
                double theta = (double)(jj + 0.5)/Nbins * M_PI;

                // Bin N particles
                AvgPhotonState avgPhoton = avgAndBinNParticles(omega, theta, 0, Nparticles, 0, 
                    everyNRGHist, finalNRGHist, everyThetaHist, finalThetaHist, everyBetaHist, finalCountHist,
                    nrgXnrgHist2D, thetaXnrgHist2D, nrgXthetaHist2D, thetaXthetaHist2D, finalValsHist2D);

                // Write these averages to output file with the current hist coordinates
                outputFile << omega << "," << theta << "," << avgPhoton.omega << "," 
                    << avgPhoton.theta << "," << avgPhoton.numScatterings << "," << avgPhoton.polarization << "\n";

                loadingBar(Nbins, ii, jj);
            }
        }

        // Close avg file
        outputFile.close();

        // Export histograms
        everyNRGHist.exportToFile("nrg", folderName);
        everyThetaHist.exportToFile("theta", folderName);
        everyBetaHist.exportToFile("beta", folderName);

        finalNRGHist.exportToFile("esc_nrg", folderName);
        finalThetaHist.exportToFile("esc_theta", folderName);
        finalCountHist.exportToFile("num", folderName);

        nrgXnrgHist2D.exportToFile("nrgXnrg", folderName);
        nrgXthetaHist2D.exportToFile("nrgXtheta", folderName);
        thetaXnrgHist2D.exportToFile("thetaXnrg", folderName);
        thetaXthetaHist2D.exportToFile("thetaXtheta", folderName);

        finalValsHist2D.exportToFile("finalVals", folderName);
    }
}