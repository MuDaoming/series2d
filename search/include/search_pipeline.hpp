#pragma once

#include <string>
#include <vector>

#include "coefficient_relation_expander.hpp"
#include "ff_type.hpp"
#include "fi_reduction_searcher.hpp"
#include "io.hpp"
#include "relation_formatter.hpp"
#include "relation_searcher.hpp"

void runRelationSearchPipeline(const std::vector<std::string>& seriesPaths,
                               const std::string& configPath,
                               const std::string& targetPath,
                               const std::string& maxSearchDegPath,
                               const std::string& outputPath);

#include "../src/search_pipeline.tpp"
