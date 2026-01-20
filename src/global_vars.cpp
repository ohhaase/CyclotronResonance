#include "global_vars.hpp"

#include <cmath>

// consts need to be set as extern as well

// ===== Simulation Parameters =====
extern const double SigmaBMin = 10;
double Theta = 0.05; // Dimensionless temp: kb * T / mc^2
extern const double wB = 0.1;
extern const int maxNum = 200; // Max number of scatters

extern const double finalNRGlow = 0.04;
extern const double finalNRGhigh = 0.2;

// ===== Constants =====
extern const double alphaf = 1.0/137.0;
extern const double B = wB;
extern const double Gamma = 4.0 * alphaf * B * B / 3.0 ;

extern const double lowerOmega = wB * sqrt((1.0 + 2.0*SigmaBMin - sqrt(1.0 + 8.0*SigmaBMin))/(2.0*(SigmaBMin - 1.0)));
extern const double upperOmega = wB * sqrt((1.0 + 2.0*SigmaBMin + sqrt(1.0 + 8.0*SigmaBMin))/(2.0*(SigmaBMin - 1.0)));

// ===== Global Variables =====
int Nbins;
int simType;
int Nparticles;
int Nthreads;

