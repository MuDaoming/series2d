#include <sstream>

template<typename T>
void RelationFormatter<T>::writeSummary(
    std::ostream& out,
    const RelationSearchResult<T>& result) {
    int pivotCount = 0;
    for (int col : result.pivotColumns) {
        if (col >= 0) {
            ++pivotCount;
        }
    }
    out << "# variables = " << result.variables.size() << "\n";
    out << "# pivot columns = " << pivotCount << "\n";
    out << "# free columns = " << result.freeColumns.size() << "\n";
}

template<typename T>
void RelationFormatter<T>::writeRelations(
    std::ostream& out,
    const RelationSearchResult<T>& result) {
    out << "[relations]\n";
    for (size_t row = 0; row < result.pivotColumns.size(); ++row) {
        const int pivotCol = result.pivotColumns[row];
        if (pivotCol < 0) {
            continue;
        }

        out << relationVariableToString(result.variables[pivotCol]) << " = ";
        bool first = true;
        for (size_t col = static_cast<size_t>(pivotCol + 1); col < result.variables.size(); ++col) {
            const T coeff = result.rrefMatrix[row][col];
            if (coeff == T(0)) {
                continue;
            }

            if (!first) {
                out << " + ";
            }
            std::ostringstream oss;
            oss << (T(0) - coeff);
            out << oss.str() << "*" << relationVariableToString(result.variables[col]);
            first = false;
        }
        if (first) {
            out << "0";
        }
        out << "\n";
    }
}

template<typename T>
void RelationFormatter<T>::writeRREF(
    std::ostream& out,
    const RelationSearchResult<T>& result) {
    out << "[rref]\n";
    for (const auto& row : result.rrefMatrix) {
        out << "{";
        for (size_t j = 0; j < row.size(); ++j) {
            std::ostringstream oss;
            oss << row[j];
            out << oss.str();
            if (j + 1 < row.size()) {
                out << ",";
            }
        }
        out << "}\n";
    }
}

template<typename T>
void RelationFormatter<T>::writeAssignments(
    std::ostream& out,
    const std::vector<CoefficientAssignment<T>>& assignments) {
    out << "[assignments]\n";
    for (const auto& assignment : assignments) {
        out << "# free column = " << assignment.chosenFreeColumn << "\n";
        for (size_t i = 0; i < assignment.variables.size(); ++i) {
            std::ostringstream oss;
            oss << assignment.values[i];
            out << relationVariableToString(assignment.variables[i])
                << " = " << oss.str() << "\n";
        }
    }
}

template<typename T>
void RelationFormatter<T>::writeFIRelations(
    std::ostream& out,
    const std::vector<FIRelation<T>>& relations) {
    out << "[fi_relations]\n";
    for (const auto& relation : relations) {
        bool first = true;
        for (size_t i = 0; i < relation.integrals.size(); ++i) {
            if (relation.coeffs[i] == T(0)) {
                continue;
            }
            if (!first) {
                out << " + ";
            }
            std::ostringstream oss;
            oss << relation.coeffs[i];
            out << oss.str() << "*" << fiVariableToString(relation.integrals[i]);
            first = false;
        }
        if (first) {
            out << "0";
        }
        out << " = 0\n";
    }
}

template<typename T>
void RelationFormatter<T>::writeFIReductionSummary(
    std::ostream& out,
    const FIReductionResult<T>& result) {
    int pivotCount = 0;
    for (int col : result.pivotColumns) {
        if (col >= 0) {
            ++pivotCount;
        }
    }
    out << "# FI variables = " << result.integrals.size() << "\n";
    out << "# FI pivot columns = " << pivotCount << "\n";
    out << "# FI free columns = " << result.freeColumns.size() << "\n";
}

template<typename T>
void RelationFormatter<T>::writeFIMasterBasis(
    std::ostream& out,
    const FIReductionResult<T>& result) {
    out << "#MIs\n";
    for (size_t i = 0; i < result.freeColumns.size(); ++i) {
        const int col = result.freeColumns[i];
        if (col < 0 || col >= static_cast<int>(result.integrals.size())) {
            continue;
        }
        out << nuToString(result.integrals[col].nu);
        if (i + 1 < result.freeColumns.size()) {
            out << ",";
        }
    }
    out << "\n";
}

template<typename T>
void RelationFormatter<T>::writeFIReductions(
    std::ostream& out,
    const FIReductionResult<T>& result) {
    out << "[fi_reductions]\n";
    for (size_t row = 0; row < result.pivotColumns.size(); ++row) {
        const int pivotCol = result.pivotColumns[row];
        if (pivotCol < 0) {
            continue;
        }

        out << fiVariableToString(result.integrals[pivotCol]) << " = ";
        bool first = true;
        for (size_t col = static_cast<size_t>(pivotCol + 1); col < result.integrals.size(); ++col) {
            const T coeff = result.rrefMatrix[row][col];
            if (coeff == T(0)) {
                continue;
            }
            if (!first) {
                out << " + ";
            }
            std::ostringstream oss;
            oss << (T(0) - coeff);
            out << oss.str() << "*" << fiVariableToString(result.integrals[col]);
            first = false;
        }
        if (first) {
            out << "0";
        }
        out << "\n";
    }
}

template<typename T>
void RelationFormatter<T>::writeFIRREF(
    std::ostream& out,
    const FIReductionResult<T>& result) {
    out << "[fi_rref]\n";
    for (const auto& row : result.rrefMatrix) {
        out << "{";
        for (size_t j = 0; j < row.size(); ++j) {
            std::ostringstream oss;
            oss << row[j];
            out << oss.str();
            if (j + 1 < row.size()) {
                out << ",";
            }
        }
        out << "}\n";
    }
}
