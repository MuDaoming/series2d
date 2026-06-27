#pragma once

#include <map>
#include <vector>

#include "absorbing_relation_searcher.hpp"
#include "relation_types.hpp"

template<typename T>
class IncrementalAbsorbingRelationSearcher {
public:
    explicit IncrementalAbsorbingRelationSearcher(const SearchInput<T>& input);
    IncrementalAbsorbingRelationSearcher(
        const SearchInput<T>& input,
        const AbsorbingSearchOptions& options);

    AbsorbingSearchResult<T> search() const;

private:
    struct PreparedData {
        int numBC = 0;
        int numIntegrals = 0;
        int degreeD = 0;
        int trainDegree = 0;
        int trainRows = 0;
        std::vector<T> coeffs;
    };

    struct RelationWork {
        std::map<int, T> coeffByVar;
    };

    struct BasisColumn {
        int pivotRow = -1;
        std::vector<T> values;
        RelationWork expression;
    };

    const SearchInput<T>& input_;
    AbsorbingSearchOptions options_;

    int varId(int integralId, int deltaPower) const;
    int varIntegral(int varId) const;
    int varPower(int varId) const;

    PreparedData prepareData() const;
    const T& coeffAt(const PreparedData& data, int bc, int integralId, int degree) const;
    std::vector<int> buildVariableOrder() const;
    std::vector<int> buildComplexOrder() const;

    std::vector<T> buildTrainColumn(const PreparedData& data, int integralId, int deltaPower) const;

    bool isAbsorbed(
        int integralId,
        int deltaPower,
        int maxM,
        const std::vector<std::vector<LeadingRule>>& rulesByIntegral) const;

    bool reduceByAccepted(
        RelationWork& relation,
        int maxM,
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

    void subtractScaled(RelationWork& dst, const RelationWork& src, const T& factor) const;
    int leadingVar(const RelationWork& relation, const std::vector<int>& orderRank) const;
    RelationWork shiftRelation(const SparsePolyRelation<T>& relation, int shift) const;
    bool reduceByRelationBasis(
        RelationWork& relation,
        const std::vector<RelationWork>& relationBasisByPivot,
        const std::vector<int>& orderRank) const;
    void insertIntoRelationBasis(
        RelationWork relation,
        std::vector<RelationWork>& relationBasisByPivot,
        const std::vector<int>& orderRank) const;
    void addRelationShiftsToBasis(
        const SparsePolyRelation<T>& relation,
        int maxM,
        std::vector<RelationWork>& relationBasisByPivot,
        const std::vector<int>& orderRank) const;
};

#include "../src/incremental_absorbing_relation_searcher.tpp"
