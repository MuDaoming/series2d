#pragma once

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <ginac/ginac.h>

struct InputConfig {
    int N = 0;
    int B = 0;
    int deg = 0;
    mp_limb_t p = 0;
    mp_limb_t a = 0;
    mp_limb_t b = 0;
    mp_limb_t d = 0;
    std::string reduceMode = "normal";
    std::string print2DMode = "target";
    std::vector<int> sector;
    std::vector<mp_limb_t> bc;
};

struct TargetConfig {
    std::vector<std::vector<int>> nus;
};

InputConfig parseConfigFile(const std::string& path);
TargetConfig parseTargetFile(const std::string& path, int expectedNuSize);
std::vector<std::vector<GiNaC::ex>> parseMatrixFile(
    const std::string& filename, GiNaC::symbol& X, GiNaC::symbol& Y);

#include "../src/io.tpp"
