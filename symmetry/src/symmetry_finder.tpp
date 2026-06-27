#pragma once

#include <algorithm>
#include <map>
#include <stdexcept>

inline SymmetryFinder::SymmetryFinder(
    const std::vector<std::vector<GiNaC::ex>>& topS,
    int numProps,
    int numBranches,
    const GiNaC::symbol& X,
    const GiNaC::symbol& Y,
    const GiNaC::ex& shiftA,
    const GiNaC::ex& shiftB)
    : numProps_(numProps),
      numBranches_(numBranches),
      X_(X),
      Y_(Y),
      shiftA_(shiftA),
      shiftB_(shiftB),
      branchOfProp_(numProps, -1),
      R_(numProps, std::vector<GiNaC::ex>(numProps)),
      canonicalizer_(R_, branchOfProp_, X, Y, shiftA, shiftB),
      transforms_(X, Y, shiftA, shiftB) {
    if (numBranches_ != 3) {
        throw std::invalid_argument("SymmetryFinder currently requires B=3");
    }
    if (topS.size() != static_cast<size_t>(numBranches_ + numProps_)) {
        throw std::invalid_argument("topS dimension does not equal B+N");
    }
    for (const auto& row : topS) {
        if (row.size() != topS.size()) {
            throw std::invalid_argument("topS must be square");
        }
    }

    for (int prop = 0; prop < numProps_; ++prop) {
        for (int branch = 0; branch < numBranches_; ++branch) {
            if (GiNaC::normal(topS[numBranches_ + prop][branch] - 1).is_zero()) {
                if (branchOfProp_[prop] != -1) {
                    throw std::runtime_error("Propagator belongs to multiple branches");
                }
                branchOfProp_[prop] = branch;
            }
        }
        if (branchOfProp_[prop] == -1) {
            throw std::runtime_error("Cannot determine propagator branch");
        }
        for (int other = 0; other < numProps_; ++other) {
            R_[prop][other] = topS[numBranches_ + prop][numBranches_ + other];
        }
    }
}

inline SectorId SymmetryFinder::sectorFromIndex(int index, int numProps) {
    SectorId result;
    result.bits.resize(numProps);
    for (int i = numProps - 1; i >= 0; --i) {
        result.bits[i] = index & 1;
        index >>= 1;
    }
    return result;
}

inline std::vector<SectorId> SymmetryFinder::enumerateValidSectors() const {
    std::vector<SectorId> result;
    const int total = 1 << numProps_;
    for (int index = total - 1; index >= 0; --index) {
        SectorId sector = sectorFromIndex(index, numProps_);
        std::vector<bool> present(numBranches_, false);
        for (int prop = 0; prop < numProps_; ++prop) {
            if (sector.bits[prop]) present[branchOfProp_[prop]] = true;
        }
        if (std::all_of(present.begin(), present.end(), [](bool value) { return value; })) {
            result.push_back(std::move(sector));
        }
    }
    return result;
}

inline SectorMapping SymmetryFinder::mappingFromWitnesses(
    const SectorCanonicalForm& source,
    const CanonicalWitness& sourceWitness,
    const SectorCanonicalForm& target,
    const CanonicalWitness& targetWitness) const {
    if (sourceWitness.orderedProps.size() != targetWitness.orderedProps.size()) {
        throw std::invalid_argument("Witness sizes differ");
    }

    SectorMapping mapping;
    mapping.source = source.sector;
    mapping.target = target.sector;
    mapping.sourceToTarget.assign(numProps_, -1);
    for (size_t i = 0; i < sourceWitness.orderedProps.size(); ++i) {
        mapping.sourceToTarget[sourceWitness.orderedProps[i]] =
            targetWitness.orderedProps[i];
    }

    for (const std::vector<int>& tau : transforms_.allBranchPermutations()) {
        bool branchCompatible = true;
        for (int sourceProp = 0; sourceProp < numProps_; ++sourceProp) {
            const int targetProp = mapping.sourceToTarget[sourceProp];
            if (targetProp < 0) continue;
            if (branchOfProp_[sourceProp] != tau[branchOfProp_[targetProp]]) {
                branchCompatible = false;
                break;
            }
        }
        if (!branchCompatible) continue;

        mapping.branchPermutation = tau;
        const auto transform = transforms_.inducedTransform(tau);
        mapping.transformedX = transform.first;
        mapping.transformedY = transform.second;
        if (verify(mapping)) return mapping;
    }
    mapping.verified = false;
    return mapping;
}

inline bool SymmetryFinder::verify(SectorMapping& mapping) const {
    if (mapping.source.bits.size() != static_cast<size_t>(numProps_) ||
        mapping.target.bits.size() != static_cast<size_t>(numProps_) ||
        mapping.sourceToTarget.size() != static_cast<size_t>(numProps_) ||
        mapping.branchPermutation.size() != static_cast<size_t>(numBranches_)) {
        mapping.verified = false;
        return false;
    }

    std::vector<int> sourceActive = activePropagators(mapping.source);
    std::vector<int> targetActive = activePropagators(mapping.target);
    if (sourceActive.size() != targetActive.size()) {
        mapping.verified = false;
        return false;
    }

    std::vector<bool> targetUsed(numProps_, false);
    for (int sourceProp : sourceActive) {
        const int targetProp = mapping.sourceToTarget[sourceProp];
        if (targetProp < 0 || targetProp >= numProps_ ||
            !mapping.target.bits[targetProp] || targetUsed[targetProp]) {
            mapping.verified = false;
            return false;
        }
        targetUsed[targetProp] = true;
        if (branchOfProp_[sourceProp] !=
            mapping.branchPermutation[branchOfProp_[targetProp]]) {
            mapping.verified = false;
            return false;
        }
    }

    const auto transform = transforms_.inducedTransform(mapping.branchPermutation);
    const GiNaC::lst substitutions = {
        X_ == transform.first,
        Y_ == transform.second
    };

    for (int sourceI : sourceActive) {
        const int targetI = mapping.sourceToTarget[sourceI];
        for (int sourceJ : sourceActive) {
            const int targetJ = mapping.sourceToTarget[sourceJ];
            const GiNaC::ex rhs = R_[targetI][targetJ].subs(
                substitutions, GiNaC::subs_options::algebraic);
            if (!GiNaC::normal(R_[sourceI][sourceJ] - rhs).is_zero()) {
                mapping.verified = false;
                return false;
            }
        }
    }

    mapping.transformedX = transform.first;
    mapping.transformedY = transform.second;
    mapping.verified = true;
    return true;
}

inline std::vector<SymmetryOrbit> SymmetryFinder::findOrbits() const {
    std::map<std::string, std::vector<SectorCanonicalForm>> groups;
    for (const SectorId& sector : enumerateValidSectors()) {
        SectorCanonicalForm form = canonicalizer_.canonicalize(sector);
        groups[form.key].push_back(std::move(form));
    }

    std::vector<SymmetryOrbit> result;
    for (auto& [key, forms] : groups) {
        std::sort(forms.begin(), forms.end(), [](const auto& lhs, const auto& rhs) {
            return sectorIndex(lhs.sector) > sectorIndex(rhs.sector);
        });

        SymmetryOrbit orbit;
        orbit.representative = forms.front().sector;
        for (const auto& form : forms) orbit.members.push_back(form.sector);

        const SectorCanonicalForm& representative = forms.front();
        for (size_t i = 1; i < forms.size(); ++i) {
            bool found = false;
            for (const CanonicalWitness& sourceWitness : forms[i].witnesses) {
                for (const CanonicalWitness& targetWitness : representative.witnesses) {
                    SectorMapping mapping = mappingFromWitnesses(
                        forms[i], sourceWitness, representative, targetWitness);
                    if (mapping.verified) {
                        orbit.mappingsToRepresentative.push_back(std::move(mapping));
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
            if (!found) {
                throw std::runtime_error(
                    "Canonical forms matched but no exact mapping verified for " +
                    sectorToString(forms[i].sector) + " -> " +
                    sectorToString(representative.sector));
            }
        }

        for (size_t i = 0; i < representative.witnesses.size(); ++i) {
            for (size_t j = i + 1; j < representative.witnesses.size(); ++j) {
                SectorMapping mapping = mappingFromWitnesses(
                    representative, representative.witnesses[i],
                    representative, representative.witnesses[j]);
                if (mapping.verified) orbit.automorphisms.push_back(std::move(mapping));
            }
        }
        result.push_back(std::move(orbit));
    }

    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        return sectorIndex(lhs.representative) > sectorIndex(rhs.representative);
    });
    return result;
}

