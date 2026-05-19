#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ginac/ginac.h>

#include "converter.hpp"
#include "delta_projector.hpp"
#include "integrand_expander.hpp"
#include "io.hpp"
#include "series_solver.hpp"

void runFI1DSeriesPipeline(const std::string& sPath,
                           const std::string& configPath,
                           const std::string& targetPath,
                           const std::string& outputPath);

#include "../src/fi_pipeline.tpp"
