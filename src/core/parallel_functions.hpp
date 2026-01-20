#pragma once

#include "../global_vars.hpp"
#include "../helper_objects/Histogram.hpp"
#include "../helper_objects/Histogram2D.hpp"

AvgPhotonState averageNParticles(double omega, double theta, int polarization, int recoil, int Nparticles);

void binNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, Histogram& totalNrgHist, Histogram& totalBetaHist);

AvgPhotonState avgAndBinNParticles(double omega, double theta, int polarization, int Nparticles, int recoil, 
    Histogram& everyNRGHist, Histogram& finalNRGHist, Histogram& everyThetaHist, Histogram& finalThetaHist,
    Histogram& everyBetaHist, Histogram& finalCountHist, Histogram2D& nrgXnrgHist2D, Histogram2D& thetaXnrgHist2D,
    Histogram2D& nrgXthetaHist2D, Histogram2D& thetaXthetaHist2D, Histogram2D& finalValsHist2D);