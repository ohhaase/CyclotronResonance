#pragma once

// ===== Cross sections =====
extern const double SigmaBMax;

double SigmaB(double w);

double cumSigmaB(double w);
extern double cumSigBMin;
extern double cumSigBMax;


// -- Polarization averaged (0) --
double averageMagTotalCrossSec(double mu_i, double w_i);
double averageMagTotalCrossSecMax(double mu_i);
double averageMagDiffCrossSec(double mu_f, double mu_i, double w_i);


// -- Perpendicular to Perpendicular (1) --
double pr2prTotalCrossSec(double mu_i, double w_i);
double pr2prTotalCrossSecMax(double mu_i);
double pr2prDiffCrossSec(double mu_f, double mu_i, double w_i);


// -- Perpendicular to Parallel (2) --
double pr2plTotalCrossSec(double mu_i, double w_i);
double pr2plTotalCrossSecMax(double mu_i);
double pr2plDiffCrossSec(double mu_f, double mu_i, double w_i);


// -- Parallel to Perpendicular (3) --
double pl2prTotalCrossSec(double mu_i, double w_i);
double pl2prTotalCrossSecMax(double mu_i);
double pl2prDiffCrossSec(double mu_f, double mu_i, double w_i);


// -- Parallel to Parallel (4) --
double pl2plTotalCrossSec(double mu_i, double w_i);
double pl2plTotalCrossSecMax(double mu_i);
double pl2plDiffCrossSec(double mu_f, double mu_i, double w_i);