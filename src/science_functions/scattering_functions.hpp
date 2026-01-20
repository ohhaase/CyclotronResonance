#pragma once

#include "../global_vars.hpp"

// ===== Scattering Functions =====
bool hasEscaped(PhotonState& photon);

double performScatter(PhotonState& photon, int recoil);