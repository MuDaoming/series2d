#pragma once

#include <algorithm>
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

struct IntegralLabel {
    IntegralHead head = IntegralHead::FI;
    std::vector<BoundaryTag> boundaries;
    std::vector<int> nu;
};

struct SeriesLabel {
    IntegralLabel integral;
    int bcIndex = 0;
};

template<typename T>
struct SeriesSample {
    SeriesLabel label;
    std::vector<T> coeffs;
};

template<typename T>
struct SearchInput {
    int degreeD = 0;
    int maxDeltaDegreeM = 0;
    int numFBIMasters = 0;
    std::vector<IntegralLabel> targets;
    std::vector<SeriesSample<T>> samples;
};

struct RelationVariable {
    IntegralLabel integral;
    int k = 0;
};

template<typename T>
struct RelationSearchResult {
    std::vector<RelationVariable> variables;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};

template<typename T>
struct CoefficientAssignment {
    std::vector<RelationVariable> variables;
    std::vector<T> values;
    int chosenFreeColumn = -1;
};

template<typename T>
struct IntegralRelation {
    std::vector<IntegralLabel> integrals;
    std::vector<T> coeffs;
};

template<typename T>
struct IntegralReductionResult {
    std::vector<IntegralLabel> integrals;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};

int countProps(const std::vector<int>& nu);
int countDots(const std::vector<int>& nu);
bool equalNu(const std::vector<int>& lhs, const std::vector<int>& rhs);
bool equalBoundaryTag(const BoundaryTag& lhs, const BoundaryTag& rhs);
bool equalIntegralLabel(const IntegralLabel& lhs, const IntegralLabel& rhs);
std::string nuToString(const std::vector<int>& nu);
std::string integralHeadToString(IntegralHead head);
std::string boundaryTagToString(const BoundaryTag& tag);
std::string integralLabelToString(const IntegralLabel& label);
std::string relationVariableToString(const RelationVariable& var);
std::string integralVariableToString(const IntegralLabel& label);

struct IntegralLabelLess {
    bool operator()(const IntegralLabel& lhs,
                    const IntegralLabel& rhs) const;
};

struct RelationVariableMoreComplexFirst {
    bool operator()(const RelationVariable& lhs,
                    const RelationVariable& rhs) const;
};

struct IntegralLabelMoreComplexFirst {
    bool operator()(const IntegralLabel& lhs,
                    const IntegralLabel& rhs) const;
};

#include "../src/relation_types.tpp"
