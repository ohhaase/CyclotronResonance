#pragma once

#include "../global_vars.hpp"
#include "../helper_objects/Histogram.hpp"
#include "../helper_objects/Histogram2D.hpp"
#include "helper_objects/OutputHandler.hpp"
#include "helper_objects/PhotonDistb.hpp"
#include "nlohmann/json.hpp"


AvgPhotonState averageNParticles(double omega, double theta, int polarization, int recoil, int Nparticles);

void binNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, Histogram& totalNrgHist, Histogram& totalBetaHist);

AvgPhotonState avgAndBinNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, 
    Histogram& everyNRGHist, Histogram& finalNRGHist, Histogram& everyThetaHist, Histogram& finalThetaHist,
    Histogram& everyBetaHist, Histogram& finalCountHist, Histogram2D& nrgXnrgHist2D, Histogram2D& thetaXnrgHist2D,
    Histogram2D& nrgXthetaHist2D, Histogram2D& thetaXthetaHist2D, Histogram2D& finalValsHist2D);


void NParticlesUniform(int Nparticles, int recoil, OutputHandler& output);

void NParticlesDistb(int Nparticles, int recoil, nlohmann::json inputDistbJSON, OutputHandler& output);