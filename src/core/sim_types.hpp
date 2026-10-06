#pragma once

#include "nlohmann/json.hpp"

// ===== Helpers =====
void storeSimInfo();
void exportSimInfo();

// ===== Sim functions =====
void simType1();
void simType2();
void simType3();
void simType4();
void simType5();
void simType6();
void simType7(nlohmann::json outputjson);
void simType8(nlohmann::json simParams);