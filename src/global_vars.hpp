#pragma once

// ===== Simulation Parameters =====
extern const double SigmaBMin;
extern double Theta; // Dimensionless temp: kb * T / mc^2
extern const double wB;
extern const int maxNum; // Max number of scatters

extern const double finalNRGlow;
extern const double finalNRGhigh;

// ===== Constants =====
extern const double alphaf;
extern const double B;
extern const double Gamma;

extern const double lowerOmega;
extern const double upperOmega;

// ===== Global Variables =====
extern int Nbins;
extern int simType;
extern int Nparticles;
extern int Nthreads;

// ===== Cross section type-defs =====
typedef double(*diffCrossSecFunc)(double, double, double); 
typedef double(*crossSecFunc)(double, double);

// ===== Data structures =====
struct PhotonState
{
    double omega;
    double theta;
    int numScatterings;
    int polarization; // 0: parallel, 1: perp
};

// Need this to have distinction between doubles and ints
struct AvgPhotonState
{
    double omega;
    double theta;
    double numScatterings;
    double polarization;
};