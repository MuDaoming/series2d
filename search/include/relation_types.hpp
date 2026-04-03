#pragma once

#include <string>
#include <vector>

struct IntegralLabel {
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
struct FIRelation {
    std::vector<IntegralLabel> integrals;
    std::vector<T> coeffs;
};

template<typename T>
struct FIReductionResult {
    std::vector<IntegralLabel> integrals;
    std::vector<std::vector<T>> rrefMatrix;
    std::vector<int> pivotColumns;
    std::vector<int> freeColumns;
};

int countProps(const std::vector<int>& nu);
int countDots(const std::vector<int>& nu);
bool equalNu(const std::vector<int>& lhs, const std::vector<int>& rhs);
std::string nuToString(const std::vector<int>& nu);
std::string relationVariableToString(const RelationVariable& var);
std::string fiVariableToString(const IntegralLabel& label);

struct RelationVariableMoreComplexFirst {
    bool operator()(const RelationVariable& lhs,
                    const RelationVariable& rhs) const;
};

struct IntegralLabelMoreComplexFirst {
    bool operator()(const IntegralLabel& lhs,
                    const IntegralLabel& rhs) const;
};

#include "../src/relation_types.tpp"
