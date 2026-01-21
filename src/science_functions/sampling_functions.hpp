#pragma once

#include "../global_vars.hpp"

// ===== Sampling functions =====
double sampleDiffCrossSection(double mu_i, double w_i, diffCrossSecFunc diffCrossSec);

double sampleSigB();

double crossSecAcceptReject(double mu_i, double targetOmega, crossSecFunc crossSec, double max);