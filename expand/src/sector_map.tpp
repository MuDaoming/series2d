#pragma once

#include <cctype>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace sector_map_detail {

inline std::string trim(const std::string& text) {
    size_t begin = 0;
    while (begin < text.size() &&
           std::isspace(static_cast<unsigned char>(text[begin]))) {
        ++begin;
    }
    size_t end = text.size();
    while (end > begin &&
           std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    return text.substr(begin, end - begin);
}

inline std::vector<int> parseList(const std::string& body) {
    std::vector<int> values;
    std::stringstream stream(body);
    std::string token;
    while (std::getline(stream, token, ',')) {
        token = trim(token);
        if (!token.empty()) values.push_back(std::stoi(token));
    }
    return values;
}

}  // namespace sector_map_detail

inline SectorMap SectorMap::fromFile(const std::string& path, int numProps) {
    std::ifstream input(path);
    if (!input.is_open()) {
        throw std::runtime_error("Cannot open sector map file: " + path);
    }

    SectorMap result(numProps);
    const std::regex linePattern(
        R"(^\s*\{\s*\{([^{}]*)\}\s*->\s*\{([^{}]*)\}\s*,\s*\{([^{}]*)\}\s*\}\s*$)");
    const std::regex pairPattern(R"((\d+)\s*->\s*(\d+))");

    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        line = sector_map_detail::trim(line);
        if (line.empty() || line[0] == '#') continue;

        std::smatch match;
        if (!std::regex_match(line, match, linePattern)) {
            throw std::runtime_error(
                "Invalid sector map line " + std::to_string(lineNumber) +
                " in " + path + ": " + line);
        }

        const std::vector<int> source =
            sector_map_detail::parseList(match[1].str());
        const std::vector<int> target =
            sector_map_detail::parseList(match[2].str());

        std::vector<std::pair<int, int>> explicitMap;
        const std::string mapBody = match[3].str();
        for (std::sregex_iterator it(mapBody.begin(), mapBody.end(), pairPattern), end;
             it != end; ++it) {
            explicitMap.emplace_back(
                std::stoi((*it)[1].str()),
                std::stoi((*it)[2].str()));
        }
        result.addEntry(source, target, explicitMap);
    }
    return result;
}

inline const SectorMapEntry* SectorMap::findEntry(
    const std::vector<int>& sector) const {
    for (const auto& entry : entries_) {
        if (entry.source == sector) return &entry;
    }
    return nullptr;
}

inline bool SectorMap::hasSource(const std::vector<int>& sector) const {
    return findEntry(sector) != nullptr;
}

inline void SectorMap::addEntry(
    const std::vector<int>& source,
    const std::vector<int>& target,
    const std::vector<std::pair<int, int>>& explicitMap) {
    if (static_cast<int>(source.size()) != numProps_ ||
        static_cast<int>(target.size()) != numProps_) {
        throw std::runtime_error("Sector map length does not match N");
    }
    if (findEntry(source)) {
        throw std::runtime_error("Duplicate source sector in sector map");
    }

    SectorMapEntry entry;
    entry.source.resize(numProps_);
    entry.target.resize(numProps_);
    entry.sourceToTarget.assign(numProps_, -1);
    for (int i = 0; i < numProps_; ++i) {
        if ((source[i] != 0 && source[i] != 1) ||
            (target[i] != 0 && target[i] != 1)) {
            throw std::runtime_error("Sector map sectors must contain only 0 or 1");
        }
        entry.source[i] = source[i];
        entry.target[i] = target[i];
    }

    std::vector<bool> targetUsed(numProps_, false);
    for (const auto& mapping : explicitMap) {
        const int sourcePos = mapping.first - 1;
        const int targetPos = mapping.second - 1;
        if (sourcePos < 0 || sourcePos >= numProps_ ||
            targetPos < 0 || targetPos >= numProps_) {
            throw std::runtime_error("Propagator position in sector map is out of range");
        }
        if (!entry.source[sourcePos] || !entry.target[targetPos]) {
            throw std::runtime_error(
                "Propagator map must connect active source and target positions");
        }
        if (entry.sourceToTarget[sourcePos] >= 0 || targetUsed[targetPos]) {
            throw std::runtime_error("Propagator map is not one-to-one");
        }
        entry.sourceToTarget[sourcePos] = targetPos;
        targetUsed[targetPos] = true;
    }

    for (int i = 0; i < numProps_; ++i) {
        if (entry.source[i] && entry.sourceToTarget[i] < 0 &&
            entry.target[i] && !targetUsed[i]) {
            entry.sourceToTarget[i] = i;
            targetUsed[i] = true;
        }
    }

    for (int sourcePos = 0; sourcePos < numProps_; ++sourcePos) {
        if (entry.source[sourcePos] && entry.sourceToTarget[sourcePos] < 0) {
            throw std::runtime_error(
                "Sector map omits a non-identity active propagator mapping");
        }
    }
    for (int targetPos = 0; targetPos < numProps_; ++targetPos) {
        if (entry.target[targetPos] && !targetUsed[targetPos]) {
            throw std::runtime_error(
                "Sector map does not cover every active target propagator");
        }
    }

    entries_.push_back(std::move(entry));
}

inline std::vector<int> SectorMap::canonicalizeSector(
    const std::vector<int>& sector) const {
    if (static_cast<int>(sector.size()) != numProps_) {
        throw std::runtime_error("Sector size does not match sector map");
    }
    std::vector<int> current = sector;
    for (size_t step = 0; step <= entries_.size(); ++step) {
        const SectorMapEntry* entry = findEntry(current);
        if (!entry) return current;
        current = entry->target;
    }
    throw std::runtime_error("Cycle detected in sector map");
}

inline std::vector<int> SectorMap::canonicalizeNu(
    const std::vector<int>& nu) const {
    if (static_cast<int>(nu.size()) != numProps_) {
        throw std::runtime_error("Nu size does not match sector map");
    }

    std::vector<int> current = nu;
    for (size_t step = 0; step <= entries_.size(); ++step) {
        std::vector<int> sector(numProps_, 0);
        for (int i = 0; i < numProps_; ++i) {
            sector[i] = current[i] > 0 ? 1 : 0;
        }
        const SectorMapEntry* entry = findEntry(sector);
        if (!entry) return current;

        std::vector<int> mapped(numProps_, 0);
        for (int sourcePos = 0; sourcePos < numProps_; ++sourcePos) {
            if (!entry->source[sourcePos]) continue;
            mapped[entry->sourceToTarget[sourcePos]] = current[sourcePos];
        }
        current = std::move(mapped);
    }
    throw std::runtime_error("Cycle detected while mapping nu");
}
