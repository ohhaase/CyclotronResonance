#pragma once

#include "../global_vars.hpp"

double calcBeta(PhotonState& photon, double w_i, double sign);

double betaForScatter(PhotonState& photon, double w_i);

bool radicalConditionMet(PhotonState& photon, double w_i);