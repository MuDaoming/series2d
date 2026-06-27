#pragma once

#include <map>
#include <vector>

#include "relation_types.hpp"

struct AbsorbVariable {
    int integralId = -1;
    int deltaPower = 0;
};

template<typename T>
struct AbsorbTerm {
    int integralId = -1;
    int deltaPower = 0;
    T coeff;
};

template<typename T>
struct SparsePolyRelation {
    std::vector<AbsorbTerm<T>> terms;
    int leadIntegralId = -1;
    int leadDeltaPower = -1;
    int maxDeltaPower = -1;
};

struct LeadingRule {
    int integralId = -1;
    int leadDeltaPower = -1;
    int relationId = -1;
    int relationMaxDeltaPower = -1;
};

struct AbsorbingSearchOptions {
    int trainDegree = -1;
    int checkStart = 0;
    int checkEnd = -1;
};

template<typename T>
struct AbsorbingSearchResult {
    std::vector<IntegralLabel> integrals;
    std::vector<SparsePolyRelation<T>> relations;
    std::vector<LeadingRule> leadingRules;
    int maxDeltaDegreeM = 0;
    int trainDegree = 0;
    int checkStart = 0;
    int checkEnd = -1;
};

template<typename T>
class AbsorbingRelationSearcher {
public:
    explicit AbsorbingRelationSearcher(const SearchInput<T>& input);
    AbsorbingRelationSearcher(
        const SearchInput<T>& input,
        const AbsorbingSearchOptions& options);

    AbsorbingSearchResult<T> search() const;

private:
    struct PreparedData {
        int numBC = 0;
        int numIntegrals = 0;
        int degreeD = 0;
        std::vector<T> coeffs;
    };

    struct RelationWork {
        std::map<int, T> coeffByVar;
    };

    const SearchInput<T>& input_;
    AbsorbingSearchOptions options_;

    int varId(int integralId, int deltaPower) const;
    int varIntegral(int varId) const;
    int varPower(int varId) const;
    std::vector<int> buildVariableOrder(int maxM) const;
    PreparedData prepareData() const;
    const T& coeffAt(const PreparedData& data, int bc, int integralId, int degree) const;

    bool isAbsorbed(
        int integralId,
        int deltaPower,
        int currentM,
        const std::vector<std::vector<LeadingRule>>& rulesByIntegral) const;

    std::vector<int> buildActiveVariables(
        int currentM,
        const std::vector<std::vector<LeadingRule>>& rulesByIntegral) const;

    std::vector<std::vector<T>> buildMatrix(
        const PreparedData& data,
        const std::vector<int>& activeVarIds,
        int trainDegree) const;

    std::vector<RelationWork> nullspaceBasis(
        const std::vector<std::vector<T>>& rref,
        const std::vector<int>& pivotColumns,
        const std::vector<int>& freeColumns,
        const std::vector<int>& activeVarIds) const;

    bool reduceByAccepted(
        RelationWork& relation,
        int currentM,
        const std::vector<SparsePolyRelation<T>>& accepted,
        const std::vector<std::vector<LeadingRule>>& rulesByIntegral,
        const std::vector<int>& orderRank) const;

    bool checkRelation(
        const RelationWork& relation,
        const PreparedData& data,
        int checkStart,
        int checkEnd) const;

    SparsePolyRelation<T> normalizeRelation(
        const RelationWork& relation,
        const std::vector<int>& orderRank) const;
};

#include "../src/absorbing_relation_searcher.tpp"
