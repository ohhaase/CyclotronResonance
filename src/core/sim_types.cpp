#include "sim_types.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <vector>

#include "nlohmann/json.hpp"

#include "../global_vars.hpp"
#include "../helper_objects/Timer.hpp"
#include "parallel_functions.hpp"
#include "helper_functions.hpp"
#include "../helper_objects/ElectronDistb.hpp"


// ===== Helpers =====
nlohmann::json simInfo;
Timer simTimer;

void storeSimInfo()
{
    simTimer.reset();

    simInfo = {
        {"StartTime", getDateTime()},
        {"RunTime", 0.0}, // Temp value, will change when we write to it
        {"NBins", Nbins},
        {"NParticles", Nparticles},
        {"NThreads", Nthreads},
        {"Recoil", true},
        {"FieldStrength", B},
        {"ElectronDistb", electronDistb.distb},
        {"ElectronTemp", electronDistb.T},
        {"NRGLims", SigmaBMin}
    };
}


void exportSimInfo()
{
    simInfo["RunTime"] = simTimer.elapsedSec();

    std::string fileName = simFolder + "/simInfo.json";

    std::ofstream simInfoFile(fileName);

    simInfoFile << std::setw(3) << simInfo << std::endl;

    simInfoFile.close();
}


// ===== Sim functions =====
void simType1()
{
    storeSimInfo();

    std::filesystem::create_directory(simFolder);
    std::ofstream outputFile(simFolder + "/" + "outdata.csv");

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

    exportSimInfo();

    std::cout << "Simulation completed in: " << totalTime.elapsedSec() << " seconds\n";

    system("py dataprocessing.py");
}

void simType2()
{
    storeSimInfo();

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
    std::filesystem::create_directory(simFolder);

    globalNRGHist.exportToFile("energy", simFolder);
    globalBetaHist.exportToFile("beta", simFolder);

    exportSimInfo();

    system("py dataprocessing.py");
}

void simType3()
{
    storeSimInfo();

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

    std::filesystem::create_directory(simFolder);
    binTimings.exportToFile("binTimings", simFolder);

    exportSimInfo();
}


void simType4()
{
    // Outputs both average and histogram data

    // Legacy version: performs 4 simulations:
    // Theta = 0.05, recoil on
    // Theta = 0.05, recoil off
    // Theta = 0.025, recoil on
    // Theta = 0.025, recoil off

    // New version: Performs as N * 2 sims, where:
    // N is number of temps we want to study
    // 2 because we look at recoil on/off for each

    std::vector<double> Thetas = {0.05, 0.025, 0.01, 0.005};
    int recoils[2] = {0, 1}; 

    // Loop over electron temeperatures
    for (int i = 0; i < Thetas.size(); i++)
    {
        electronDistb.updateTemp(Thetas[i]);

        // Loop over recoil
        for (int j = 0; j < 2; j++)
        {
            storeSimInfo();
            simInfo["Recoil"] = (bool)recoils[j];

            // Unique folder
            simFolder = "simData" + std::to_string(j + i*2);
            std::filesystem::create_directory(simFolder);

            // Initialize output file for average data
            std::ofstream outputFile(simFolder + "/" + "avgdata.csv");

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
            everyNRGHist.exportToFile("nrg", simFolder);
            everyThetaHist.exportToFile("theta", simFolder);
            everyBetaHist.exportToFile("beta", simFolder);

            finalNRGHist.exportToFile("esc_nrg", simFolder);
            finalThetaHist.exportToFile("esc_theta", simFolder);
            finalCountHist.exportToFile("num", simFolder);

            nrgXnrgHist2D.exportToFile("nrgXnrg", simFolder);
            nrgXthetaHist2D.exportToFile("nrgXtheta", simFolder);
            thetaXnrgHist2D.exportToFile("thetaXnrg", simFolder);
            thetaXthetaHist2D.exportToFile("thetaXtheta", simFolder);

            finalValsHist2D.exportToFile("finalVals", simFolder);

            exportSimInfo();
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
        electronDistb.updateTemp(Thetas[i]);

        // Loop over recoil
        for (int j = 0; j < 2; j++)
        {
            storeSimInfo();
            simInfo["Recoil"] = (bool)recoils[j];

            // Unique folder
            simFolder = "simData" + std::to_string(j + i*2);
            std::filesystem::create_directory(simFolder);

            // Initialize output file for average data
            std::ofstream outputFile(simFolder + "/" + "avgdata.csv");

            // Initialize histograms
            Histogram everyNRGHist{Nbins, lowerOmega, upperOmega};
            Histogram finalNRGHist{Nbins, finalNRGlow, finalNRGhigh, true};
            Histogram everyThetaHist{Nbins, 0, M_PI};
            Histogram finalThetaHist{Nbins, 0, M_PI};
            Histogram everyBetaHist{Nbins, -1, 1};
            Histogram finalCountHist{Nbins, 1, maxNum};

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
            everyNRGHist.exportToFile("nrg", simFolder);
            everyThetaHist.exportToFile("theta", simFolder);
            everyBetaHist.exportToFile("beta", simFolder);

            finalNRGHist.exportToFile("esc_nrg", simFolder);
            finalThetaHist.exportToFile("esc_theta", simFolder);
            finalCountHist.exportToFile("num", simFolder);

            nrgXnrgHist2D.exportToFile("nrgXnrg", simFolder);
            nrgXthetaHist2D.exportToFile("nrgXtheta", simFolder);
            thetaXnrgHist2D.exportToFile("thetaXnrg", simFolder);
            thetaXthetaHist2D.exportToFile("thetaXtheta", simFolder);

            finalValsHist2D.exportToFile("finalVals", simFolder);

            exportSimInfo();
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
        electronDistb.updateTemp(Thetas[i]);
        storeSimInfo();

        // Unique folder
        simFolder = "simData" + std::to_string(i);
        std::filesystem::create_directory(simFolder);

        // Initialize output file for average data
        std::ofstream outputFile(simFolder + "/" + "avgdata.csv");

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
        everyNRGHist.exportToFile("nrg", simFolder);
        everyThetaHist.exportToFile("theta", simFolder);
        everyBetaHist.exportToFile("beta", simFolder);

        finalNRGHist.exportToFile("esc_nrg", simFolder);
        finalThetaHist.exportToFile("esc_theta", simFolder);
        finalCountHist.exportToFile("num", simFolder);

        nrgXnrgHist2D.exportToFile("nrgXnrg", simFolder);
        nrgXthetaHist2D.exportToFile("nrgXtheta", simFolder);
        thetaXnrgHist2D.exportToFile("thetaXnrg", simFolder);
        thetaXthetaHist2D.exportToFile("thetaXtheta", simFolder);

        finalValsHist2D.exportToFile("finalVals", simFolder);

        exportSimInfo();
    }
}