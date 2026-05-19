#pragma once

#include <stdexcept>
#include <string>
#include <vector>

enum class IntegralHead {
    FI,
    BFI,
    BBFI
};

enum class BoundaryAxis {
    X,
    Y
};

enum class BoundarySide {
    U,
    D
};

struct BoundaryTag {
    BoundaryAxis axis;
    BoundarySide side;
};

struct IntegralTag {
    IntegralHead head = IntegralHead::FI;
    std::vector<BoundaryTag> boundaries;
    std::vector<int> nu;
};

inline std::string toString(IntegralHead head) {
    switch (head) {
        case IntegralHead::FI: return "FI";
        case IntegralHead::BFI: return "BFI";
        case IntegralHead::BBFI: return "BBFI";
    }
    throw std::runtime_error("Unknown IntegralHead");
}

inline std::string toString(const BoundaryTag& tag) {
    std::string s;
    s += (tag.axis == BoundaryAxis::X) ? "X" : "Y";
    s += (tag.side == BoundarySide::U) ? "U" : "D";
    return s;
}
