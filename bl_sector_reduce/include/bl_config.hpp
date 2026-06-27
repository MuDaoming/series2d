#pragma once

#include <flint/nmod.h>

struct BLSectorConfig {
    int nuSize = 0;
    int degreeD = 0;
    int maxDegree = 0;
    int safetyOrder = 10;
    int certOrder = 10;
    mp_limb_t prime = 0;
};

