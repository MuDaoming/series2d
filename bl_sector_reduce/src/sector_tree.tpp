#pragma once

#include <algorithm>
#include <cctype>
#include <queue>
#include <sstream>
#include <stdexcept>

bool equalSectorId(const SectorId& lhs, const SectorId& rhs) {
    return lhs.bits == rhs.bits;
}

bool sectorContains(const SectorId& lhs, const SectorId& rhs) {
    if (lhs.bits.size() != rhs.bits.size()) return false;
    for (size_t i = 0; i < lhs.bits.size(); ++i) {
        if (rhs.bits[i] && !lhs.bits[i]) return false;
    }
    return true;
}

int sectorPopcount(const SectorId& sector) {
    int n = 0;
    for (int b : sector.bits) if (b) ++n;
    return n;
}

std::string sectorIdToString(const SectorId& sector) {
    std::string s = "{";
    for (size_t i = 0; i < sector.bits.size(); ++i) {
        s += std::to_string(sector.bits[i]);
        if (i + 1 < sector.bits.size()) s += ",";
    }
    s += "}";
    return s;
}

SectorId parseSectorId(const std::string& text, int expectedNuSize) {
    const size_t l = text.find('{');
    const size_t r = text.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Invalid sector id: " + text);
    }
    SectorId id;
    std::stringstream ss(text.substr(l + 1, r - l - 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        size_t a = 0;
        while (a < tok.size() && std::isspace(static_cast<unsigned char>(tok[a]))) ++a;
        size_t b = tok.size();
        while (b > a && std::isspace(static_cast<unsigned char>(tok[b - 1]))) --b;
        if (a == b) continue;
        int bit = std::stoi(tok.substr(a, b - a));
        id.bits.push_back(bit ? 1 : 0);
    }
    if (expectedNuSize >= 0 && static_cast<int>(id.bits.size()) != expectedNuSize) {
        throw std::runtime_error("Sector length mismatch: " + text);
    }
    return id;
}

bool SectorIdLess::operator()(const SectorId& lhs, const SectorId& rhs) const {
    return std::lexicographical_compare(
        lhs.bits.begin(), lhs.bits.end(), rhs.bits.begin(), rhs.bits.end());
}

SectorTree::SectorTree(const std::vector<SectorId>& sectors) : sectors_(sectors) {
    std::sort(sectors_.begin(), sectors_.end(), SectorIdLess{});
    sectors_.erase(std::unique(sectors_.begin(), sectors_.end(), equalSectorId), sectors_.end());
    if (sectors_.empty()) throw std::runtime_error("SectorTree needs at least one sector");
    buildOrder();
}

int SectorTree::rootIndex() const { return rootIndex_; }
const std::vector<int>& SectorTree::processingOrder() const { return order_; }
int SectorTree::parentOf(int sectorIndex) const { return parents_.at(sectorIndex); }
const SectorId& SectorTree::sectorAt(int idx) const { return sectors_.at(idx); }
int SectorTree::size() const { return static_cast<int>(sectors_.size()); }

int SectorTree::indexOf(const SectorId& sector) const {
    for (int i = 0; i < static_cast<int>(sectors_.size()); ++i) {
        if (equalSectorId(sectors_[i], sector)) return i;
    }
    throw std::runtime_error("Unknown sector: " + sectorIdToString(sector));
}

std::vector<int> SectorTree::ancestorsOf(int sectorIndex) const {
    std::vector<int> out;
    const SectorId& sector = sectors_.at(sectorIndex);
    const int sectorComplexity = sectorPopcount(sector);
    for (int idx : order_) {
        if (idx == sectorIndex) break;
        if (sectorPopcount(sectors_[idx]) > sectorComplexity) {
            out.push_back(idx);
        }
    }
    return out;
}

void SectorTree::buildOrder() {
    const int n = static_cast<int>(sectors_.size());
    parents_.assign(n, -1);
    order_.resize(n);
    for (int i = 0; i < n; ++i) order_[i] = i;
    std::sort(order_.begin(), order_.end(), [&](int a, int b) {
        const int pa = sectorPopcount(sectors_[a]);
        const int pb = sectorPopcount(sectors_[b]);
        if (pa != pb) return pa > pb;
        return SectorIdLess{}(sectors_[a], sectors_[b]);
    });
    rootIndex_ = order_.front();
}
