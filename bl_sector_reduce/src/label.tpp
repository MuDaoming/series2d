#pragma once

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

static std::string blTrim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

static ObjectHead parseObjectHead(const std::string& raw) {
    if (raw == "FI") return ObjectHead::FI;
    if (raw == "BFI") return ObjectHead::BFI;
    if (raw == "BBFI") return ObjectHead::BBFI;
    throw std::runtime_error("Invalid object head: " + raw);
}

static BoundaryTag parseBoundaryTag(const std::string& raw) {
    const std::string s = blTrim(raw);
    if (s.size() != 2) {
        throw std::runtime_error("Invalid boundary tag: " + raw);
    }
    BoundaryTag tag;
    if (s[0] == 'X') tag.axis = BoundaryAxis::X;
    else if (s[0] == 'Y') tag.axis = BoundaryAxis::Y;
    else throw std::runtime_error("Invalid boundary axis: " + raw);

    if (s[1] == 'U') tag.side = BoundarySide::U;
    else if (s[1] == 'D') tag.side = BoundarySide::D;
    else throw std::runtime_error("Invalid boundary side: " + raw);
    return tag;
}

static std::vector<BoundaryTag> parseBoundaryList(const std::string& raw) {
    std::vector<BoundaryTag> tags;
    std::stringstream ss(raw);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = blTrim(tok);
        if (!tok.empty()) tags.push_back(parseBoundaryTag(tok));
    }
    return tags;
}

static void validateObjectLabel(const ObjectLabel& label, const std::string& line) {
    if (label.head == ObjectHead::FI) {
        if (!label.boundaries.empty()) {
            throw std::runtime_error("FI object must not have boundaries: " + line);
        }
        return;
    }
    if (label.head == ObjectHead::BFI) {
        if (label.boundaries.size() != 1) {
            throw std::runtime_error("BFI object must have exactly one boundary: " + line);
        }
        return;
    }
    if (label.boundaries.size() != 2) {
        throw std::runtime_error("BBFI object must have exactly two boundaries: " + line);
    }
    bool hasX = false;
    bool hasY = false;
    for (const auto& b : label.boundaries) {
        hasX = hasX || b.axis == BoundaryAxis::X;
        hasY = hasY || b.axis == BoundaryAxis::Y;
    }
    if (!hasX || !hasY) {
        throw std::runtime_error("BBFI object must have one X and one Y boundary: " + line);
    }
}

bool equalBoundaryTag(const BoundaryTag& lhs, const BoundaryTag& rhs) {
    return lhs.axis == rhs.axis && lhs.side == rhs.side;
}

bool equalObjectLabel(const ObjectLabel& lhs, const ObjectLabel& rhs) {
    if (lhs.head != rhs.head || lhs.nu != rhs.nu ||
        lhs.boundaries.size() != rhs.boundaries.size()) {
        return false;
    }
    for (size_t i = 0; i < lhs.boundaries.size(); ++i) {
        if (!equalBoundaryTag(lhs.boundaries[i], rhs.boundaries[i])) return false;
    }
    return true;
}

ObjectLabel parseObjectLabel(const std::string& text, int expectedNuSize) {
    const std::string line = blTrim(text);
    const size_t l = line.find('{');
    const size_t r = line.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Invalid object label: " + text);
    }

    ObjectLabel label;
    const std::string prefix = blTrim(line.substr(0, l));
    if (prefix.empty()) {
        label.head = ObjectHead::FI;
    } else {
        const size_t lb = prefix.find('[');
        if (lb == std::string::npos) {
            label.head = parseObjectHead(prefix);
        } else {
            const size_t rb = prefix.rfind(']');
            if (rb == std::string::npos || rb <= lb || rb + 1 != prefix.size()) {
                throw std::runtime_error("Invalid boundary list in label: " + text);
            }
            label.head = parseObjectHead(blTrim(prefix.substr(0, lb)));
            label.boundaries = parseBoundaryList(prefix.substr(lb + 1, rb - lb - 1));
        }
    }

    std::stringstream ss(line.substr(l + 1, r - l - 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = blTrim(tok);
        if (!tok.empty()) label.nu.push_back(std::stoi(tok));
    }
    if (expectedNuSize >= 0 && static_cast<int>(label.nu.size()) != expectedNuSize) {
        throw std::runtime_error("nu size mismatch in object label: " + text);
    }
    validateObjectLabel(label, line);
    return label;
}

std::string objectHeadToString(ObjectHead head) {
    switch (head) {
        case ObjectHead::FI: return "FI";
        case ObjectHead::BFI: return "BFI";
        case ObjectHead::BBFI: return "BBFI";
    }
    throw std::runtime_error("Unknown ObjectHead");
}

std::string boundaryTagToString(const BoundaryTag& tag) {
    std::string s;
    s += tag.axis == BoundaryAxis::X ? "X" : "Y";
    s += tag.side == BoundarySide::U ? "U" : "D";
    return s;
}

std::string objectLabelToString(const ObjectLabel& label) {
    std::string s = objectHeadToString(label.head);
    if (!label.boundaries.empty()) {
        s += "[";
        for (size_t i = 0; i < label.boundaries.size(); ++i) {
            s += boundaryTagToString(label.boundaries[i]);
            if (i + 1 < label.boundaries.size()) s += ",";
        }
        s += "]";
    }
    s += "{";
    for (size_t i = 0; i < label.nu.size(); ++i) {
        s += std::to_string(label.nu[i]);
        if (i + 1 < label.nu.size()) s += ",";
    }
    s += "}";
    return s;
}

std::vector<int> sectorSupport(const ObjectLabel& label) {
    std::vector<int> bits;
    bits.reserve(label.nu.size());
    for (int x : label.nu) bits.push_back(x > 0 ? 1 : 0);
    return bits;
}

static int boundaryRank(const BoundaryTag& tag) {
    const int axis = tag.axis == BoundaryAxis::X ? 0 : 1;
    const int side = tag.side == BoundarySide::U ? 0 : 1;
    return 2 * axis + side;
}

bool ObjectLabelLess::operator()(const ObjectLabel& lhs, const ObjectLabel& rhs) const {
    if (lhs.head != rhs.head) {
        return static_cast<int>(lhs.head) < static_cast<int>(rhs.head);
    }
    if (lhs.boundaries.size() != rhs.boundaries.size()) {
        return lhs.boundaries.size() < rhs.boundaries.size();
    }
    for (size_t i = 0; i < lhs.boundaries.size(); ++i) {
        const int lr = boundaryRank(lhs.boundaries[i]);
        const int rr = boundaryRank(rhs.boundaries[i]);
        if (lr != rr) return lr < rr;
    }
    return std::lexicographical_compare(
        lhs.nu.begin(), lhs.nu.end(), rhs.nu.begin(), rhs.nu.end());
}
