#include <algorithm>
#include <iostream>
#include <stdexcept>

template<typename T>
IncrementalAbsorbingRelationSearcher<T>::IncrementalAbsorbingRelationSearcher(
    const SearchInput<T>& input)
    : input_(input) {
    options_.trainDegree = input.degreeD;
    options_.checkStart = 0;
    options_.checkEnd = input.degreeD;
}

template<typename T>
IncrementalAbsorbingRelationSearcher<T>::IncrementalAbsorbingRelationSearcher(
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
        throw std::runtime_error("incremental absorbing search degree window exceeds input degree");
    }
}

template<typename T>
int IncrementalAbsorbingRelationSearcher<T>::varId(int integralId, int deltaPower) const {
    return integralId * (input_.maxDeltaDegreeM + 1) + deltaPower;
}

template<typename T>
int IncrementalAbsorbingRelationSearcher<T>::varIntegral(int id) const {
    return id / (input_.maxDeltaDegreeM + 1);
}

template<typename T>
int IncrementalAbsorbingRelationSearcher<T>::varPower(int id) const {
    return id % (input_.maxDeltaDegreeM + 1);
}

template<typename T>
typename IncrementalAbsorbingRelationSearcher<T>::PreparedData
IncrementalAbsorbingRelationSearcher<T>::prepareData() const {
    if (input_.numFBIMasters <= 0) {
        throw std::runtime_error("numFBIMasters must be positive");
    }
    const int nInt = static_cast<int>(input_.targets.size());
    if (nInt <= 0) {
        throw std::runtime_error("incremental absorbing search requires nonempty targets");
    }
    if (static_cast<int>(input_.samples.size()) != input_.numFBIMasters * nInt) {
        throw std::runtime_error("sample count must be numFBIMasters * target count");
    }

    PreparedData data;
    data.numBC = input_.numFBIMasters;
    data.numIntegrals = nInt;
    data.degreeD = input_.degreeD;
    data.trainDegree = options_.trainDegree;
    data.trainRows = data.numBC * (data.trainDegree + 1);
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
const T& IncrementalAbsorbingRelationSearcher<T>::coeffAt(
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
std::vector<int> IncrementalAbsorbingRelationSearcher<T>::buildVariableOrder() const {
    std::vector<int> ids;
    ids.reserve(input_.targets.size() * static_cast<size_t>(input_.maxDeltaDegreeM + 1));
    for (int i = 0; i < static_cast<int>(input_.targets.size()); ++i) {
        for (int k = 0; k <= input_.maxDeltaDegreeM; ++k) {
            ids.push_back(varId(i, k));
        }
    }
    std::sort(ids.begin(), ids.end(), [&](int a, int b) {
        RelationVariable va{input_.targets[varIntegral(a)], varPower(a)};
        RelationVariable vb{input_.targets[varIntegral(b)], varPower(b)};
        RelationVariableMoreComplexFirst moreComplex;
        if (moreComplex(va, vb)) return false;
        if (moreComplex(vb, va)) return true;
        return a < b;
    });
    return ids;
}

template<typename T>
std::vector<int> IncrementalAbsorbingRelationSearcher<T>::buildComplexOrder() const {
    std::vector<int> ids;
    ids.reserve(input_.targets.size() * static_cast<size_t>(input_.maxDeltaDegreeM + 1));
    for (int i = 0; i < static_cast<int>(input_.targets.size()); ++i) {
        for (int k = 0; k <= input_.maxDeltaDegreeM; ++k) {
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
std::vector<T> IncrementalAbsorbingRelationSearcher<T>::buildTrainColumn(
    const PreparedData& data,
    int integralId,
    int deltaPower) const {
    std::vector<T> col(static_cast<size_t>(data.trainRows), T(0));
    int row = 0;
    for (int bc = 0; bc < data.numBC; ++bc) {
        for (int n = 0; n <= data.trainDegree; ++n) {
            if (n >= deltaPower) {
                col[static_cast<size_t>(row)] =
                    coeffAt(data, bc, integralId, n - deltaPower);
            }
            ++row;
        }
    }
    return col;
}

template<typename T>
bool IncrementalAbsorbingRelationSearcher<T>::isAbsorbed(
    int integralId,
    int deltaPower,
    int maxM,
    const std::vector<std::vector<LeadingRule>>& rulesByIntegral) const {
    for (const auto& rule : rulesByIntegral[static_cast<size_t>(integralId)]) {
        if (deltaPower < rule.leadDeltaPower) {
            continue;
        }
        const int shift = deltaPower - rule.leadDeltaPower;
        if (rule.relationMaxDeltaPower + shift <= maxM) {
            return true;
        }
    }
    return false;
}

template<typename T>
void IncrementalAbsorbingRelationSearcher<T>::subtractScaled(
    RelationWork& dst,
    const RelationWork& src,
    const T& factor) const {
    const T zero(0);
    for (const auto& kv : src.coeffByVar) {
        T& out = dst.coeffByVar[kv.first];
        out = out - factor * kv.second;
        if (out == zero) {
            dst.coeffByVar.erase(kv.first);
        }
    }
}

template<typename T>
int IncrementalAbsorbingRelationSearcher<T>::leadingVar(
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
    return leadVar;
}

template<typename T>
typename IncrementalAbsorbingRelationSearcher<T>::RelationWork
IncrementalAbsorbingRelationSearcher<T>::shiftRelation(
    const SparsePolyRelation<T>& relation,
    int shift) const {
    RelationWork out;
    for (const auto& term : relation.terms) {
        out.coeffByVar[varId(term.integralId, term.deltaPower + shift)] = term.coeff;
    }
    return out;
}

template<typename T>
bool IncrementalAbsorbingRelationSearcher<T>::reduceByRelationBasis(
    RelationWork& relation,
    const std::vector<RelationWork>& relationBasisByPivot,
    const std::vector<int>& orderRank) const {
    while (!relation.coeffByVar.empty()) {
        const int lead = leadingVar(relation, orderRank);
        if (lead < 0) {
            relation.coeffByVar.clear();
            return false;
        }
        const auto& reducer = relationBasisByPivot[static_cast<size_t>(lead)];
        if (reducer.coeffByVar.empty()) {
            return true;
        }
        const T factor = relation.coeffByVar[lead];
        subtractScaled(relation, reducer, factor);
    }
    return false;
}

template<typename T>
void IncrementalAbsorbingRelationSearcher<T>::insertIntoRelationBasis(
    RelationWork relation,
    std::vector<RelationWork>& relationBasisByPivot,
    const std::vector<int>& orderRank) const {
    if (!reduceByRelationBasis(relation, relationBasisByPivot, orderRank)) {
        return;
    }
    const int lead = leadingVar(relation, orderRank);
    if (lead < 0) {
        return;
    }
    const T invLead = T(1) / relation.coeffByVar[lead];
    for (auto& kv : relation.coeffByVar) {
        kv.second = kv.second * invLead;
    }
    relationBasisByPivot[static_cast<size_t>(lead)] = std::move(relation);
}

template<typename T>
void IncrementalAbsorbingRelationSearcher<T>::addRelationShiftsToBasis(
    const SparsePolyRelation<T>& relation,
    int maxM,
    std::vector<RelationWork>& relationBasisByPivot,
    const std::vector<int>& orderRank) const {
    for (int shift = 0; relation.maxDeltaPower + shift <= maxM; ++shift) {
        insertIntoRelationBasis(
            shiftRelation(relation, shift),
            relationBasisByPivot,
            orderRank);
    }
}

template<typename T>
bool IncrementalAbsorbingRelationSearcher<T>::reduceByAccepted(
    RelationWork& relation,
    int maxM,
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
            if (rule.relationMaxDeltaPower + shift <= maxM) {
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
bool IncrementalAbsorbingRelationSearcher<T>::checkRelation(
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
SparsePolyRelation<T> IncrementalAbsorbingRelationSearcher<T>::normalizeRelation(
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
AbsorbingSearchResult<T> IncrementalAbsorbingRelationSearcher<T>::search() const {
    const PreparedData data = prepareData();
    const int maxM = input_.maxDeltaDegreeM;

    const auto variableOrder = buildVariableOrder();
    const auto complexOrder = buildComplexOrder();
    std::vector<int> orderRank(input_.targets.size() * static_cast<size_t>(maxM + 1), 0);
    for (int rank = 0; rank < static_cast<int>(complexOrder.size()); ++rank) {
        orderRank[static_cast<size_t>(complexOrder[static_cast<size_t>(rank)])] = rank;
    }

    std::vector<BasisColumn> basis;
    std::vector<int> pivotToBasis(static_cast<size_t>(data.trainRows), -1);
    std::vector<SparsePolyRelation<T>> accepted;
    std::vector<LeadingRule> leadingRules;
    std::vector<std::vector<LeadingRule>> rulesByIntegral(input_.targets.size());
    std::vector<RelationWork> relationBasisByPivot(
        input_.targets.size() * static_cast<size_t>(maxM + 1));

    const T zero(0);
    for (int var : variableOrder) {
        const int integralId = varIntegral(var);
        const int power = varPower(var);
        if (isAbsorbed(integralId, power, maxM, rulesByIntegral)) {
            continue;
        }

        std::vector<T> col = buildTrainColumn(data, integralId, power);
        RelationWork expr;
        expr.coeffByVar[var] = T(1);

        for (const auto& b : basis) {
            const T factor = col[static_cast<size_t>(b.pivotRow)];
            if (factor == zero) {
                continue;
            }
            for (int row = b.pivotRow; row < data.trainRows; ++row) {
                col[static_cast<size_t>(row)] =
                    col[static_cast<size_t>(row)] -
                    factor * b.values[static_cast<size_t>(row)];
            }
            subtractScaled(expr, b.expression, factor);
        }

        int pivot = -1;
        for (int row = 0; row < data.trainRows; ++row) {
            if (col[static_cast<size_t>(row)] != zero) {
                pivot = row;
                break;
            }
        }

        if (pivot >= 0) {
            const T invPivot = T(1) / col[static_cast<size_t>(pivot)];
            for (int row = pivot; row < data.trainRows; ++row) {
                if (col[static_cast<size_t>(row)] != zero) {
                    col[static_cast<size_t>(row)] =
                        col[static_cast<size_t>(row)] * invPivot;
                }
            }
            for (auto& kv : expr.coeffByVar) {
                kv.second = kv.second * invPivot;
            }
            BasisColumn out;
            out.pivotRow = pivot;
            out.values = std::move(col);
            out.expression = std::move(expr);
            pivotToBasis[static_cast<size_t>(pivot)] = static_cast<int>(basis.size());
            basis.push_back(std::move(out));
            continue;
        }

        if (!reduceByAccepted(expr, maxM, accepted, rulesByIntegral, orderRank)) {
            continue;
        }
        if (!reduceByRelationBasis(expr, relationBasisByPivot, orderRank)) {
            continue;
        }
        if (!checkRelation(expr, data, options_.checkStart, options_.checkEnd)) {
            continue;
        }
        auto rel = normalizeRelation(expr, orderRank);
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
        addRelationShiftsToBasis(accepted.back(), maxM, relationBasisByPivot, orderRank);
    }

    AbsorbingSearchResult<T> result;
    result.integrals = input_.targets;
    result.relations = std::move(accepted);
    result.leadingRules = std::move(leadingRules);
    result.maxDeltaDegreeM = maxM;
    result.trainDegree = options_.trainDegree;
    result.checkStart = options_.checkStart;
    result.checkEnd = options_.checkEnd;
    return result;
}
