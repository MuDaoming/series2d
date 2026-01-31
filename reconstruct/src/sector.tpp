// Perform RREF and compute rowOperation matrix simultaneously
template<typename T>
void Sector<T>::rowReduce() {
    int n = S_.size();
    
    // Initialize reducedS_ as copy of S_
    reducedS_ = S_;
    
    // Initialize rowOperation as identity matrix
    rowOperation_.assign(n, std::vector<T>(n, T(0)));
    for (int i = 0; i < n; ++i) {
        rowOperation_[i][i] = T(1);
    }
    
    int rank = 0;
    int current_row = 0;
    
    // Perform RREF with simultaneous operations on rowOperation
    for (int col = 0; col < n && current_row < n; ++col) {
        // Find pivot in current column
        int pivot_row = -1;
        for (int row = current_row; row < n; ++row) {
            if (reducedS_[row][col] != T(0)) {
                pivot_row = row;
                break;
            }
        }
        
        // Skip if column is all zeros
        if (pivot_row == -1) {
            continue;
        }

        rank ++;
        
        // Swap rows if needed
        if (pivot_row != current_row) {
            std::swap(reducedS_[pivot_row], reducedS_[current_row]);
            std::swap(rowOperation_[pivot_row], rowOperation_[current_row]);
        }
        
        T pivot = reducedS_[current_row][col];
        
        // Scale pivot row to make pivot = 1
        for (int j = 0; j < n; ++j) {
            reducedS_[current_row][j] = reducedS_[current_row][j] / pivot;
            rowOperation_[current_row][j] = rowOperation_[current_row][j] / pivot;
        }
        
        // Eliminate other rows (both above and below)
        for (int row = 0; row < n; ++row) {
            if (row != current_row && reducedS_[row][col] != T(0)) {
                T factor = reducedS_[row][col];
                for (int j = 0; j < n; ++j) {
                    reducedS_[row][j] = reducedS_[row][j] - factor * reducedS_[current_row][j];
                    rowOperation_[row][j] = rowOperation_[row][j] - factor * rowOperation_[current_row][j];
                }
            }
        }
        
        current_row++;
    }
    
    // Compute dimNull_
    // int rank = 0;
    // for (int row = 0; row < n; ++row) {
    //     bool is_zero_row = true;
    //     for (int col = 0; col < n; ++col) {
    //         if (reducedS_[row][col] != T(0)) {
    //             is_zero_row = false;
    //             break;
    //         }
    //     }
    //     if (!is_zero_row) {
    //         rank++;
    //     }
    // }
    dimNull_ = n - rank;
    
    // Set z0_
    z0_ = (dimNull_ == 0) ? T(1) : T(0);
}

// Find null space basis (for dimNull_ >= 1)
template<typename T>
std::vector<std::vector<T>> Sector<T>::findNullSpace() const {
    if (dimNull_ == 0) {
        return std::vector<std::vector<T>>(); // empty basis
    }
    
    int n = S_.size();
    std::vector<std::vector<T>> null_basis;
    
    // Find free variables (columns without pivots in reducedS)
    std::vector<bool> is_pivot_col(n, false);
    std::vector<int> free_vars;
    
    // Mark pivot columns
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            if (reducedS_[row][col] != T(0)) {
                is_pivot_col[col] = true;
                break; // found pivot for this row
            }
        }
    }
    
    // Collect free variables
    for (int col = 0; col < n; ++col) {
        if (!is_pivot_col[col]) {
            free_vars.push_back(col);
        }
    }
    
    // Generate basis vectors for null space
    for (int free_var : free_vars) {
        std::vector<T> null_vector(n, T(0));
        null_vector[free_var] = T(1);
        
        // Back-substitute to find dependent variables
        for (int row = n - 1; row >= 0; --row) {
            // Find pivot column for this row
            int pivot_col = -1;
            for (int col = 0; col < n; ++col) {
                if (reducedS_[row][col] != T(0)) {
                    pivot_col = col;
                    break;
                }
            }
            
            if (pivot_col != -1 && pivot_col != free_var) {
                T sum = T(0);
                for (int col = pivot_col + 1; col < n; ++col) {
                    sum += reducedS_[row][col] * null_vector[col];
                }
                null_vector[pivot_col] = -sum;
            }
        }
        
        null_basis.push_back(null_vector);
    }
    
    return null_basis;
}

// Solve linear system (for dimNull_ == 0)
template<typename T>
std::vector<T> Sector<T>::solveLinear(const std::vector<T>& b) const {
    if (dimNull_ != 0) {
        throw std::runtime_error("Matrix is not invertible (dimNull_ != 0)");
    }
    
    // For dimNull_ == 0, rowOperation is the inverse matrix
    // Solve S * x = b => x = invS * b = rowOperation * b
    int n = S_.size();
    std::vector<T> solution(n, T(0));
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n && j < b.size(); ++j) {
            solution[i] += rowOperation_[i][j] * b[j];
        }
    }
    
    return solution;
}

// Solve for c and z coefficients
template<typename T>
void Sector<T>::solveCandZ() {
    int n = S_.size(); // total matrix size (N+B)x(N+B)
    
    if (dimNull_ == 0) {
        // Case: dimNull_ == 0
        // S (C_1,...,C_B,z_1,...,z_N) = (1_1,...,1_B,0,...,0)
        std::vector<T> rhs(n, T(0));
        for (int i = 0; i < numBranch_ && i < n; ++i) {
            rhs[i] = T(1);
        }
        candz_ = solveLinear(rhs);
    } else if (dimNull_ == 1) {
        // Case: dimNull_ == 1
        std::vector<std::vector<T>> null_basis = findNullSpace();
        if (!null_basis.empty()) {
            candz_ = null_basis[0];
        }
    } else {
        // Case: dimNull_ >= 2
        std::vector<std::vector<T>> null_basis = findNullSpace();
        T best_C_sum = T(0);
        std::vector<T> best_solution(n, T(0));
        bool found_nonzero = false;
        for (const auto& basis_vector : null_basis) {
            T candidate_C_sum = T(0);
            for (int i = 0; i < numBranch_ && i < basis_vector.size(); ++i) {
                candidate_C_sum += basis_vector[i];
            }
            if (!found_nonzero || candidate_C_sum != T(0)) {
                best_C_sum = candidate_C_sum;
                best_solution = basis_vector;
                if (candidate_C_sum != T(0)) {
                    found_nonzero = true;
                    break;
                }
            }
        }
        candz_ = best_solution;
    }
    
    // Compute C_ = sum of c_i
    C_ = T(0);
    for (int i = 0; i < numBranch_; ++i) {
        C_ += candz_[i];
    }
}

