#include <algorithm>
#include <stdexcept>

#include "linear.hpp"

template<typename T>
AbsorbingRelationSearcher<T>::AbsorbingRelationSearcher(const SearchInput<T>& input)
    : input_(input) {
    options_.trainDegree = input.degreeD;
    options_.checkStart = 0;
    options_.checkEnd = input.degreeD;
}

template<typename T>
AbsorbingRelationSearcher<T>::AbsorbingRelationSearcher(
    const SearchInput<T>& input,
    const AbsorbingSearchOptions& options)
    : input_(input), options_(options) {
    if (options_.trainDegree < 0) {
        options_.trainDegree = input.degreeD;
    }
    if (options_.checkEnd < options_.checkStart) {
        options_.checkStart = 0;
        options_.checkEnd = input.degreeD;
    }
    if (options_.trainDegree > input.degreeD || options_.checkEnd > input.degreeD) {
        throw std::runtime_error("absorbing search degree window exceeds input degree");
    }
}

template<typename T>
int AbsorbingRelationSearcher<T>::varId(int integralId, int deltaPower) const {
    return integralId * (input_.maxDeltaDegreeM + 1) + deltaPower;
}

template<typename T>
int AbsorbingRelationSearcher<T>::varIntegral(int id) const {
    return id / (input_.maxDeltaDegreeM + 1);
}

template<typename T>
int AbsorbingRelationSearcher<T>::varPower(int id) const {
    return id % (input_.maxDeltaDegreeM + 1);
}

template<typename T>
std::vector<int> AbsorbingRelationSearcher<T>::buildVariableOrder(int maxM) const {
    std::vector<int> ids;
    ids.reserve(input_.targets.size() * static_cast<size_t>(maxM + 1));
    for (int i = 0; i < static_cast<int>(input_.targets.size()); ++i) {
        for (int k = 0; k <= maxM; ++k) {
            ids.push_back(varId(i, k));
        }
    }
    std::sort(ids.begin(), ids.end(), [&](int a, int b) {
        RelationVariable va{input_.targets[varIntegral(a)], varPower(a)};
        RelationVariable vb{input_.targets[varIntegral(b)], varPower(b)};
        return RelationVariableMoreComplexFirst{}(va, vb);
    });
    return ids;
}

template<typename T>
typename AbsorbingRelationSearcher<T>::PreparedData
AbsorbingRelationSearcher<T>::prepareData() const {
    if (input_.numFBIMasters <= 0) {
        throw std::runtime_error("numFBIMasters must be positive");
    }
    const int nInt = static_cast<int>(input_.targets.size());
    if (nInt <= 0) {
        throw std::runtime_error("absorbing search requires nonempty targets");
    }
    if (static_cast<int>(input_.samples.size()) != input_.numFBIMasters * nInt) {
        throw std::runtime_error("sample count must be numFBIMasters * target count");
    }

    PreparedData data;
    data.numBC = input_.numFBIMasters;
    data.numIntegrals = nInt;
    data.degreeD = input_.degreeD;
    data.coeffs.assign(
        static_cast<size_t>(data.numBC) * nInt * (data.degreeD + 1), T(0));

    for (const auto& sample : input_.samples) {
        int integralId = -1;
        for (int i = 0; i < nInt; ++i) {
            if (equalIntegralLabel(sample.label.integral, input_.targets[i])) {
                integralId = i;
                break;
            }
        }
        if (integralId < 0) {
            throw std::runtime_error("sample integral is not in targets");
        }
        if (sample.label.bcIndex < 0 || sample.label.bcIndex >= data.numBC) {
            throw std::runtime_error("sample bcIndex out of range");
        }
        if (static_cast<int>(sample.coeffs.size()) < data.degreeD + 1) {
            throw std::runtime_error("sample coeff count smaller than degreeD+1");
        }
        for (int d = 0; d <= data.degreeD; ++d) {
            data.coeffs[
                (static_cast<size_t>(sample.label.bcIndex) * nInt + integralId) *
                    (data.degreeD + 1) +
                d] = sample.coeffs[d];
        }
    }
    return data;
}

template<typename T>
const T& AbsorbingRelationSearcher<T>::coeffAt(
    const PreparedData& data,
    int bc,
    int integralId,
    int degree) const {
    return data.coeffs[
        (static_cast<size_t>(bc) * data.numIntegrals + integralId) *
            (data.degreeD + 1) +
        degree];
}

template<typename T>
bool AbsorbingRelationSearcher<T>::isAbsorbed(
    int integralId,
    int deltaPower,
    int currentM,
    const std::vector<std::vector<LeadingRule>>& rulesByIntegral) const {
    for (const auto& rule : rulesByIntegral[static_cast<size_t>(integralId)]) {
        if (deltaPower < rule.leadDeltaPower) {
            continue;
        }
        const int shift = deltaPower - rule.leadDeltaPower;
        if (rule.relationMaxDeltaPower + shift <= currentM) {
            return true;
        }
    }
    return false;
}

template<typename T>
std::vector<int> AbsorbingRelationSearcher<T>::buildActiveVariables(
    int currentM,
    const std::vector<std::vector<LeadingRule>>& rulesByIntegral) const {
    std::vector<int> active;
    active.reserve(input_.targets.size() * static_cast<size_t>(currentM + 1));
    const auto ordered = buildVariableOrder(currentM);
    for (int id : ordered) {
        const int integralId = varIntegral(id);
        const int k = varPower(id);
        if (!isAbsorbed(integralId, k, currentM, rulesByIntegral)) {
            active.push_back(id);
        }
    }
    return active;
}

template<typename T>
std::vector<std::vector<T>> AbsorbingRelationSearcher<T>::buildMatrix(
    const PreparedData& data,
    const std::vector<int>& activeVarIds,
    int trainDegree) const {
    std::vector<std::vector<T>> matrix;
    if (activeVarIds.empty()) {
        return matrix;
    }
    matrix.reserve(static_cast<size_t>(data.numBC) * (trainDegree + 1));
    for (int bc = 0; bc < data.numBC; ++bc) {
        for (int n = 0; n <= trainDegree; ++n) {
            std::vector<T> row(activeVarIds.size(), T(0));
            for (int c = 0; c < static_cast<int>(activeVarIds.size()); ++c) {
                const int id = activeVarIds[static_cast<size_t>(c)];
                const int k = varPower(id);
                if (n >= k) {
                    row[static_cast<size_t>(c)] =
                        coeffAt(data, bc, varIntegral(id), n - k);
                }
            }
            matrix.push_back(std::move(row));
        }
    }
    return matrix;
}

template<typename T>
std::vector<typename AbsorbingRelationSearcher<T>::RelationWork>
AbsorbingRelationSearcher<T>::nullspaceBasis(
    const std::vector<std::vector<T>>& rref,
    const std::vector<int>& pivotColumns,
    const std::vector<int>& freeColumns,
    const std::vector<int>& activeVarIds) const {
    std::vector<RelationWork> basis;
    basis.reserve(freeColumns.size());
    const T zero(0);
    for (int freeCol : freeColumns) {
        RelationWork rel;
        rel.coeffByVar[activeVarIds[static_cast<size_t>(freeCol)]] = T(1);
        for (int row = 0; row < static_cast<int>(pivotColumns.size()); ++row) {
            const int pivot = pivotColumns[static_cast<size_t>(row)];
            if (pivot < 0) {
                continue;
            }
            const T coeff = T(0) - rref[static_cast<size_t>(row)][static_cast<size_t>(freeCol)];
            if (coeff != zero) {
                rel.coeffByVar[activeVarIds[static_cast<size_t>(pivot)]] = coeff;
            }
        }
        basis.push_back(std::move(rel));
    }
    return basis;
}

template<typename T>
bool AbsorbingRelationSearcher<T>::reduceByAccepted(
    RelationWork& relation,
    int currentM,
    const std::vector<SparsePolyRelation<T>>& accepted,
    const std::vector<std::vector<LeadingRule>>& rulesByIntegral,
    const std::vector<int>& orderRank) const {
    const T zero(0);
    while (!relation.coeffByVar.empty()) {
        int leadVar = -1;
        int leadRank = static_cast<int>(orderRank.size()) + 1;
        for (const auto& kv : relation.coeffByVar) {
            if (kv.second == zero) {
                continue;
            }
            const int rank = orderRank[static_cast<size_t>(kv.first)];
            if (rank < leadRank) {
                leadRank = rank;
                leadVar = kv.first;
            }
        }
        if (leadVar < 0) {
            relation.coeffByVar.clear();
            return false;
        }

        const int integralId = varIntegral(leadVar);
        const int power = varPower(leadVar);
        const LeadingRule* best = nullptr;
        for (const auto& rule : rulesByIntegral[static_cast<size_t>(integralId)]) {
            if (power < rule.leadDeltaPower) {
                continue;
            }
            const int shift = power - rule.leadDeltaPower;
            if (rule.relationMaxDeltaPower + shift <= currentM) {
                best = &rule;
                break;
            }
        }
        if (best == nullptr) {
            return true;
        }

        const auto& base = accepted[static_cast<size_t>(best->relationId)];
        const int shift = power - best->leadDeltaPower;
        const T leadCoeff = relation.coeffByVar[leadVar];
        for (const auto& term : base.terms) {
            const int shiftedPower = term.deltaPower + shift;
            const int shiftedVar = varId(term.integralId, shiftedPower);
            T& dst = relation.coeffByVar[shiftedVar];
            dst = dst - leadCoeff * term.coeff;
            if (dst == zero) {
                relation.coeffByVar.erase(shiftedVar);
            }
        }
    }
    return false;
}

template<typename T>
bool AbsorbingRelationSearcher<T>::checkRelation(
    const RelationWork& relation,
    const PreparedData& data,
    int checkStart,
    int checkEnd) const {
    if (checkStart > checkEnd) {
        return true;
    }
    const T zero(0);
    for (int bc = 0; bc < data.numBC; ++bc) {
        for (int n = checkStart; n <= checkEnd; ++n) {
            T sum(0);
            for (const auto& kv : relation.coeffByVar) {
                const int id = kv.first;
                const int k = varPower(id);
                if (n >= k) {
                    sum += kv.second * coeffAt(data, bc, varIntegral(id), n - k);
                }
            }
            if (sum != zero) {
                return false;
            }
        }
    }
    return true;
}

template<typename T>
SparsePolyRelation<T> AbsorbingRelationSearcher<T>::normalizeRelation(
    const RelationWork& relation,
    const std::vector<int>& orderRank) const {
    const T zero(0);
    int leadVar = -1;
    int leadRank = static_cast<int>(orderRank.size()) + 1;
    for (const auto& kv : relation.coeffByVar) {
        if (kv.second == zero) {
            continue;
        }
        const int rank = orderRank[static_cast<size_t>(kv.first)];
        if (rank < leadRank) {
            leadRank = rank;
            leadVar = kv.first;
        }
    }
    if (leadVar < 0) {
        throw std::runtime_error("cannot normalize zero relation");
    }
    const T invLead = T(1) / relation.coeffByVar.at(leadVar);

    SparsePolyRelation<T> out;
    out.leadIntegralId = varIntegral(leadVar);
    out.leadDeltaPower = varPower(leadVar);
    out.maxDeltaPower = 0;
    out.terms.reserve(relation.coeffByVar.size());
    for (const auto& kv : relation.coeffByVar) {
        const T coeff = kv.second * invLead;
        if (coeff == zero) {
            continue;
        }
        AbsorbTerm<T> term;
        term.integralId = varIntegral(kv.first);
        term.deltaPower = varPower(kv.first);
        term.coeff = coeff;
        out.maxDeltaPower = std::max(out.maxDeltaPower, term.deltaPower);
        out.terms.push_back(term);
    }
    std::sort(out.terms.begin(), out.terms.end(), [&](const auto& a, const auto& b) {
        const int va = varId(a.integralId, a.deltaPower);
        const int vb = varId(b.integralId, b.deltaPower);
        return orderRank[static_cast<size_t>(va)] < orderRank[static_cast<size_t>(vb)];
    });
    return out;
}

template<typename T>
AbsorbingSearchResult<T> AbsorbingRelationSearcher<T>::search() const {
    const PreparedData data = prepareData();
    const int maxM = input_.maxDeltaDegreeM;
    const int trainDegree = options_.trainDegree;

    const auto fullOrder = buildVariableOrder(maxM);
    std::vector<int> orderRank(input_.targets.size() * static_cast<size_t>(maxM + 1), 0);
    for (int rank = 0; rank < static_cast<int>(fullOrder.size()); ++rank) {
        orderRank[static_cast<size_t>(fullOrder[static_cast<size_t>(rank)])] = rank;
    }

    std::vector<SparsePolyRelation<T>> accepted;
    std::vector<LeadingRule> leadingRules;
    std::vector<std::vector<LeadingRule>> rulesByIntegral(input_.targets.size());

    for (int m = 0; m <= maxM; ++m) {
        const auto activeVars = buildActiveVariables(m, rulesByIntegral);
        if (activeVars.empty()) {
            continue;
        }
        auto matrix = buildMatrix(data, activeVars, trainDegree);
        if (matrix.empty()) {
            continue;
        }
        LinearSystem<T> system(matrix);
        system.eliminate();
        auto candidates = nullspaceBasis(
            system.getRREFMatrix(),
            system.getPivotColumns(),
            system.getFreeVariableColumns(),
            activeVars);

        for (auto& candidate : candidates) {
            if (!reduceByAccepted(candidate, m, accepted, rulesByIntegral, orderRank)) {
                continue;
            }
            if (!checkRelation(candidate, data, options_.checkStart, options_.checkEnd)) {
                continue;
            }
            auto rel = normalizeRelation(candidate, orderRank);
            LeadingRule rule;
            rule.integralId = rel.leadIntegralId;
            rule.leadDeltaPower = rel.leadDeltaPower;
            rule.relationId = static_cast<int>(accepted.size());
            rule.relationMaxDeltaPower = rel.maxDeltaPower;
            accepted.push_back(std::move(rel));
            leadingRules.push_back(rule);
            auto& bucket = rulesByIntegral[static_cast<size_t>(rule.integralId)];
            bucket.push_back(rule);
            std::sort(bucket.begin(), bucket.end(), [](const auto& a, const auto& b) {
                if (a.leadDeltaPower != b.leadDeltaPower) {
                    return a.leadDeltaPower < b.leadDeltaPower;
                }
                return a.relationMaxDeltaPower < b.relationMaxDeltaPower;
            });
        }
    }

    AbsorbingSearchResult<T> result;
    result.integrals = input_.targets;
    result.relations = std::move(accepted);
    result.leadingRules = std::move(leadingRules);
    result.maxDeltaDegreeM = maxM;
    result.trainDegree = trainDegree;
    result.checkStart = options_.checkStart;
    result.checkEnd = options_.checkEnd;
    return result;
}
