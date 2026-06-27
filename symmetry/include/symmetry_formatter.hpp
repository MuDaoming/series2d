#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "symmetry_types.hpp"

void writeSymmetryReport(std::ostream& out,
                         const std::vector<SymmetryOrbit>& orbits);

#include "../src/symmetry_formatter.tpp"

