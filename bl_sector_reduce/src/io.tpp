#pragma once

#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

static std::string ioTrim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

static std::string resolveRelative(const std::string& baseFile, const std::string& path) {
    std::filesystem::path p(path);
    if (p.is_absolute()) return p.string();
    return (std::filesystem::path(baseFile).parent_path() / p).lexically_normal().string();
}

BLSectorConfig parseBLSectorConfig(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open config file: " + path);
    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = ioTrim(line);
        if (line.empty() || line[0] == '#') continue;
        const size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        kv[ioTrim(line.substr(0, eq))] = ioTrim(line.substr(eq + 1));
    }
    auto need = [&](const std::string& key) -> std::string {
        auto it = kv.find(key);
        if (it == kv.end()) throw std::runtime_error("Missing config key: " + key);
        return it->second;
    };
    BLSectorConfig cfg;
    cfg.nuSize = std::stoi(need("N"));
    cfg.degreeD = std::stoi(need("deg"));
    cfg.maxDegree = std::stoi(need("m"));
    cfg.prime = static_cast<mp_limb_t>(std::stoull(need("p")));
    if (kv.count("K_safety")) cfg.safetyOrder = std::stoi(kv["K_safety"]);
    if (kv.count("K_cert")) cfg.certOrder = std::stoi(kv["K_cert"]);
    if (cfg.safetyOrder < 0 || cfg.certOrder < 0) {
        throw std::runtime_error("K_safety and K_cert must be nonnegative");
    }
    return cfg;
}

std::vector<SectorSeriesEntry> parseSectorSeriesList(const std::string& path,
                                                     int expectedNuSize) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open sector series list: " + path);
    std::vector<SectorSeriesEntry> entries;
    std::string line;
    while (std::getline(in, line)) {
        line = ioTrim(line);
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string sectorTok;
        SectorSeriesEntry e;
        if (!(ss >> sectorTok >> e.seriesPath >> e.targetPath >> e.masterPath)) {
            throw std::runtime_error("Invalid sector series list line: " + line);
        }
        const size_t eq = sectorTok.find('=');
        if (eq == std::string::npos || sectorTok.substr(0, eq) != "sector") {
            throw std::runtime_error("Sector list line must start with sector={...}: " + line);
        }
        e.sector = parseSectorId(sectorTok.substr(eq + 1), expectedNuSize);
        e.seriesPath = resolveRelative(path, e.seriesPath);
        e.targetPath = resolveRelative(path, e.targetPath);
        e.masterPath = resolveRelative(path, e.masterPath);
        entries.push_back(std::move(e));
    }
    if (entries.empty()) throw std::runtime_error("Sector series list is empty: " + path);
    return entries;
}

std::vector<ObjectLabel> parseObjectList(const std::string& path, int expectedNuSize) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open object list: " + path);
    std::vector<ObjectLabel> out;
    std::string line;
    while (std::getline(in, line)) {
        line = ioTrim(line);
        if (line.empty() || line[0] == '#') continue;
        out.push_back(parseObjectLabel(line, expectedNuSize));
    }
    if (out.empty()) throw std::runtime_error("Object list is empty: " + path);
    return out;
}

std::vector<ObjectLabel> parseTargetFile(const std::string& path, int expectedNuSize) {
    return parseObjectList(path, expectedNuSize);
}

template<typename T>
static std::vector<T> parseSeriesLine(const std::string& line) {
    const size_t l = line.find('{');
    const size_t r = line.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Invalid series line: " + line);
    }
    std::vector<T> coeffs;
    std::stringstream ss(line.substr(l + 1, r - l - 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = ioTrim(tok);
        if (!tok.empty()) coeffs.emplace_back(static_cast<unsigned long long>(std::stoull(tok)));
    }
    return coeffs;
}

template<typename T>
void loadSectorData(const std::vector<SectorSeriesEntry>& entries,
                    int degreeD,
                    int expectedNuSize,
                    SeriesStore<T>& series,
                    MasterData& masters) {
    for (const auto& e : entries) {
        const auto targets = parseTargetFile(e.targetPath, expectedNuSize);
        std::ifstream sin(e.seriesPath);
        if (!sin.is_open()) throw std::runtime_error("Cannot open series file: " + e.seriesPath);
        std::string line;
        int idx = 0;
        while (std::getline(sin, line)) {
            line = ioTrim(line);
            if (line.empty()) continue;
            if (idx >= static_cast<int>(targets.size())) {
                throw std::runtime_error("Too many series lines in " + e.seriesPath);
            }
            auto coeffs = parseSeriesLine<T>(line);
            if (static_cast<int>(coeffs.size()) != degreeD + 1) {
                throw std::runtime_error("Series degree mismatch in " + e.seriesPath);
            }
            series.addSeries(e.sector, targets[idx], std::move(coeffs));
            ++idx;
        }
        if (idx != static_cast<int>(targets.size())) {
            throw std::runtime_error("Series/target line mismatch in " + e.seriesPath);
        }

        std::ifstream min(e.masterPath);
        if (!min.is_open()) throw std::runtime_error("Cannot open master file: " + e.masterPath);
        std::vector<ObjectLabel> ms;
        while (std::getline(min, line)) {
            line = ioTrim(line);
            if (line.empty() || line[0] == '#') continue;
            ms.push_back(parseObjectLabel(line, expectedNuSize));
        }
        masters.setMasters(e.sector, std::move(ms));
    }
}

