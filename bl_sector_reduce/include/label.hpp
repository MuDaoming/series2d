#pragma once

#include <string>
#include <vector>

enum class ObjectHead { FI, BFI, BBFI };
enum class BoundaryAxis { X, Y };
enum class BoundarySide { U, D };

struct BoundaryTag {
    BoundaryAxis axis;
    BoundarySide side;
};

struct ObjectLabel {
    ObjectHead head = ObjectHead::FI;
    std::vector<BoundaryTag> boundaries;
    std::vector<int> nu;
};

bool equalBoundaryTag(const BoundaryTag& lhs, const BoundaryTag& rhs);
bool equalObjectLabel(const ObjectLabel& lhs, const ObjectLabel& rhs);
ObjectLabel parseObjectLabel(const std::string& text, int expectedNuSize);
std::string objectHeadToString(ObjectHead head);
std::string boundaryTagToString(const BoundaryTag& tag);
std::string objectLabelToString(const ObjectLabel& label);
std::vector<int> sectorSupport(const ObjectLabel& label);

struct ObjectLabelLess {
    bool operator()(const ObjectLabel& lhs, const ObjectLabel& rhs) const;
};

#include "../src/label.tpp"
