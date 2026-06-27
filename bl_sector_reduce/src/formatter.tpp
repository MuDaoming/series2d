#pragma once

#include <ostream>
#include <sstream>

template<typename T>
std::string polynomialToString(const Polynomial1D<T>& poly) {
    std::ostringstream oss;
    oss << "{";
    const int d = std::max(0, poly.degree());
    for (int i = 0; i <= d; ++i) {
        if (i) oss << ",";
        oss << poly.coeff(i);
    }
    oss << "}";
    return oss.str();
}

template<typename T>
void writeReductions(std::ostream& out,
                     const BLSectorConfig& config,
                     const SectorTree& tree,
                     const MasterData& masters,
                     const std::vector<SectorReduction<T>>& reductions) {
    out << "# p = " << config.prime << "\n";
    out << "# D = " << config.degreeD << "\n";
    out << "# m_user = " << config.maxDegree << "\n";
    out << "# K_safety = " << config.safetyOrder << "\n";
    out << "# K_cert = " << config.certOrder << "\n\n";
    out << "# sectors = " << tree.size() << "\n";
    int masterCount = 0;
    for (int idx : tree.processingOrder()) {
        masterCount += static_cast<int>(masters.mastersFor(tree.sectorAt(idx)).size());
    }
    out << "# input masters = " << masterCount << "\n";
    out << "# reductions = " << reductions.size() << "\n\n";

    out << "[sector_reductions]\n";
    for (const auto& red : reductions) {
        out << "sector=" << sectorIdToString(red.sector) << "\n";
        out << "object=" << objectLabelToString(red.object) << "\n";
        if (red.failed) {
            out << "status=failed\n";
            out << "error=" << red.failureReason << "\n";
            if (!red.residualPrefix.empty()) {
                out << "residual_prefix={";
                for (size_t i = 0; i < red.residualPrefix.size(); ++i) {
                    if (i) out << ",";
                    out << red.residualPrefix[i];
                }
                out << "}\n";
            }
            out << "\n";
            continue;
        }
        out << "status=success\n";
        if (red.isZero) {
            out << "zero\n\n";
            continue;
        }
        out << "free_master=" << (red.isFreeMaster ? 1 : 0) << "\n";
        out << "den=" << polynomialToString(red.denominator) << "\n";
        for (const auto& term : red.terms) {
            out << "term " << objectLabelToString(term.master)
                << "=" << polynomialToString(term.numerator) << "\n";
        }
        out << "\n";
    }

    out << "[global_reductions]\n";
    for (const auto& red : reductions) {
        if (red.failed || red.isZero) continue;
        out << objectLabelToString(red.object) << " @ "
            << sectorIdToString(red.sector) << " : ";
        out << "den " << polynomialToString(red.denominator);
        for (const auto& term : red.terms) {
            out << " + (" << polynomialToString(term.numerator) << ")/den*"
                << objectLabelToString(term.master);
        }
        out << "\n";
    }
}
