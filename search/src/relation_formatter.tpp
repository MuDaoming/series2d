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
void RelationFormatter<T>::writeIntegralRelations(
    std::ostream& out,
    const std::vector<IntegralRelation<T>>& relations) {
    out << "[relations]\n";
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
            out << oss.str() << "*" << integralVariableToString(relation.integrals[i]);
            first = false;
        }
        if (first) {
            out << "0";
        }
        out << " = 0\n";
    }
}

template<typename T>
void RelationFormatter<T>::writeIntegralReductionSummary(
    std::ostream& out,
    const IntegralReductionResult<T>& result) {
    int pivotCount = 0;
    for (int col : result.pivotColumns) {
        if (col >= 0) {
            ++pivotCount;
        }
    }
    out << "# integral variables = " << result.integrals.size() << "\n";
    out << "# integral pivot columns = " << pivotCount << "\n";
    out << "# integral free columns = " << result.freeColumns.size() << "\n";
}

template<typename T>
void RelationFormatter<T>::writeIntegralMasterBasis(
    std::ostream& out,
    const IntegralReductionResult<T>& result) {
    out << "#MIs\n";
    for (size_t i = 0; i < result.freeColumns.size(); ++i) {
        const int col = result.freeColumns[i];
        if (col < 0 || col >= static_cast<int>(result.integrals.size())) {
            continue;
        }
        out << integralLabelToString(result.integrals[col]);
        if (i + 1 < result.freeColumns.size()) {
            out << ",";
        }
    }
    out << "\n";
}

template<typename T>
void RelationFormatter<T>::writeIntegralReductions(
    std::ostream& out,
    const IntegralReductionResult<T>& result) {
    out << "[reductions]\n";
    for (size_t row = 0; row < result.pivotColumns.size(); ++row) {
        const int pivotCol = result.pivotColumns[row];
        if (pivotCol < 0) {
            continue;
        }

        out << integralVariableToString(result.integrals[pivotCol]) << " = ";
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
            out << oss.str() << "*" << integralVariableToString(result.integrals[col]);
            first = false;
        }
        if (first) {
            out << "0";
        }
        out << "\n";
    }
}

template<typename T>
void RelationFormatter<T>::writeIntegralRREF(
    std::ostream& out,
    const IntegralReductionResult<T>& result) {
    out << "[integral_rref]\n";
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
