#pragma once

#include <flint/nmod.h>
#include <string>

struct BLSectorConfig {
    int nuSize = 0;
    int degreeD = -1;
    int maxDegree = 0;
    int safetyOrder = 10;
    int certOrder = 10;
    mp_limb_t prime = 0;
    std::string sectorMapPath;
};
