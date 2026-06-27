#pragma once

#include <algorithm>
#include <sstream>
#include <stdexcept>

inline SectorCanonicalizer::SectorCanonicalizer(
    const std::vector<std::vector<GiNaC::ex>>& R,
    const std::vector<int>& branchOfProp,
    const GiNaC::symbol& X,
    const GiNaC::symbol& Y,
    const GiNaC::ex& shiftA,
    const GiNaC::ex& shiftB)
    : R_(R),
      branchOfProp_(branchOfProp),
      X_(X),
      Y_(Y),
      transforms_(X, Y, shiftA, shiftB) {
    if (R_.size() != branchOfProp_.size()) {
        throw std::invalid_argument("R size and branchOfProp size differ");
    }
    for (const auto& row : R_) {
        if (row.size() != R_.size()) {
            throw std::invalid_argument("R must be square");
        }
    }
}

inline std::string SectorCanonicalizer::expressionKey(
    const GiNaC::ex& value) {
    const GiNaC::ex normalized = GiNaC::normal(value);
    const GiNaC::ex numerator = GiNaC::expand(GiNaC::numer(normalized));
    const GiNaC::ex denominator = GiNaC::expand(GiNaC::denom(normalized));
    std::ostringstream out;
    out << "(" << numerator << ")/(" << denominator << ")";
    return out.str();
}

inline std::string SectorCanonicalizer::encodeTokens(
    const std::vector<std::string>& tokens) {
    std::ostringstream out;
    for (const std::string& token : tokens) {
        out << token.size() << ":" << token;
    }
    return out.str();
}

inline std::vector<std::vector<GiNaC::ex>>
SectorCanonicalizer::transformedR(
    const std::pair<GiNaC::ex, GiNaC::ex>& transform) const {
    std::vector<std::vector<GiNaC::ex>> result = R_;
    const GiNaC::lst substitutions = {
        X_ == transform.first,
        Y_ == transform.second
    };
    for (auto& row : result) {
        for (GiNaC::ex& entry : row) {
            entry = GiNaC::normal(entry.subs(substitutions, GiNaC::subs_options::algebraic));
        }
    }
    return result;
}

inline std::vector<CanonicalWitness> SectorCanonicalizer::canonicalizeForTau(
    const SectorId& sector,
    const std::vector<int>& tau) const {
    const auto transform = transforms_.inducedTransform(tau);
    const auto matrix = transformedR(transform);
    const std::vector<int> active = activePropagators(sector);

    std::vector<int> canonicalBranches;
    for (int destinationBranch = 0; destinationBranch < 3; ++destinationBranch) {
        const int sourceBranch = tau[destinationBranch];
        for (int prop : active) {
            if (branchOfProp_[prop] == sourceBranch) {
                canonicalBranches.push_back(destinationBranch);
            }
        }
    }
    if (canonicalBranches.size() != active.size()) {
        throw std::runtime_error("Failed to build canonical branch positions");
    }

    PartialCandidate initial;
    initial.used.assign(R_.size(), false);
    std::vector<PartialCandidate> candidates = {initial};

    for (size_t pos = 0; pos < canonicalBranches.size(); ++pos) {
        const int sourceBranch = tau[canonicalBranches[pos]];
        std::vector<PartialCandidate> expanded;

        for (const PartialCandidate& candidate : candidates) {
            for (int prop : active) {
                if (candidate.used[prop] || branchOfProp_[prop] != sourceBranch) {
                    continue;
                }
                PartialCandidate next = candidate;
                for (int previous : next.orderedProps) {
                    next.prefixTokens.push_back(expressionKey(matrix[prop][previous]));
                }
                next.prefixTokens.push_back(expressionKey(matrix[prop][prop]));
                next.orderedProps.push_back(prop);
                next.used[prop] = true;
                expanded.push_back(std::move(next));
            }
        }

        if (expanded.empty()) {
            throw std::runtime_error("No branch-compatible propagator candidate");
        }

        std::string best = encodeTokens(expanded.front().prefixTokens);
        for (size_t i = 1; i < expanded.size(); ++i) {
            best = std::min(best, encodeTokens(expanded[i].prefixTokens));
        }

        candidates.clear();
        for (PartialCandidate& candidate : expanded) {
            if (encodeTokens(candidate.prefixTokens) == best) {
                candidates.push_back(std::move(candidate));
            }
        }
    }

    std::vector<CanonicalWitness> witnesses;
    for (const PartialCandidate& candidate : candidates) {
        CanonicalWitness witness;
        witness.branchPermutation = tau;
        witness.orderedProps = candidate.orderedProps;
        witness.key = encodeTokens(candidate.prefixTokens);
        witnesses.push_back(std::move(witness));
    }
    return witnesses;
}

inline SectorCanonicalForm SectorCanonicalizer::canonicalize(
    const SectorId& sector) const {
    if (sector.bits.size() != R_.size()) {
        throw std::invalid_argument("Sector size does not match R");
    }

    SectorCanonicalForm result;
    result.sector = sector;

    for (const std::vector<int>& tau : transforms_.allBranchPermutations()) {
        std::vector<CanonicalWitness> witnesses =
            canonicalizeForTau(sector, tau);
        for (CanonicalWitness& witness : witnesses) {
            if (result.witnesses.empty() || witness.key < result.key) {
                result.key = witness.key;
                result.witnesses = {std::move(witness)};
            } else if (witness.key == result.key) {
                result.witnesses.push_back(std::move(witness));
            }
        }
    }
    return result;
}
