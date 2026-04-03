#pragma once

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ff_type.hpp"
#include "relation_types.hpp"

struct SearchConfig {
    int nuSize = 0;
    int degreeD = 0;
    int numFBIMasters = 0;
    mp_limb_t p = 0;
};

SearchConfig parseSearchConfigFile(const std::string& path);
int parseMaxSearchDegreeFile(const std::string& path);
std::vector<IntegralLabel> parseSearchTargetFile(const std::string& path, int expectedNuSize);

template<typename T>
std::vector<SeriesSample<T>> parseSeriesFile(
    const std::string& path,
    const std::vector<IntegralLabel>& targets,
    int degreeD,
    int bcIndex);

#include "../src/io.tpp"
