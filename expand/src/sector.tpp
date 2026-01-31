/**
 * Sector class 实现
 * 处理子矩阵S，进行RREF行化简，计算C和z系数
 */

// Perform RREF and compute rowOperation matrix simultaneously
template<typename RT, typename PT>
void Sector<RT, PT>::rowReduce() {
    int n = S_.size();
    
    // Initialize reducedS_ as copy of S_ (PT -> RT conversion)
    reducedS_.resize(n, std::vector<RT>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            reducedS_[i][j] = RT(S_[i][j]);  // PT -> RT 转换
        }
    }
    
    // Initialize rowOperation as identity matrix
    rowOperation_.assign(n, std::vector<RT>(n, RT(0)));
    for (int i = 0; i < n; ++i) {
        rowOperation_[i][i] = RT(1);
    }
    
    int rank = 0;
    int current_row = 0;
    
    // Perform RREF with simultaneous operations on rowOperation
    for (int col = 0; col < n && current_row < n; ++col) {
        // Find pivot in current column
        int pivot_row = -1;
        for (int row = current_row; row < n; ++row) {
            if (!isZero(reducedS_[row][col])) {
                pivot_row = row;
                break;
            }
        }
        
        // Skip if column is all zeros
        if (pivot_row == -1) {
            continue;
        }

        rank++;
        
        // Swap rows if needed
        if (pivot_row != current_row) {
            std::swap(reducedS_[pivot_row], reducedS_[current_row]);
            std::swap(rowOperation_[pivot_row], rowOperation_[current_row]);
        }
        
        RT pivot = reducedS_[current_row][col];
        
        // Scale pivot row to make pivot = 1
        for (int j = 0; j < n; ++j) {
            reducedS_[current_row][j] = normalize(reducedS_[current_row][j] / pivot);
            rowOperation_[current_row][j] = normalize(rowOperation_[current_row][j] / pivot);
        }
        
        // Eliminate other rows (both above and below)
        for (int row = 0; row < n; ++row) {
            if (row != current_row && !isZero(reducedS_[row][col])) {
                RT factor = reducedS_[row][col];
                for (int j = 0; j < n; ++j) {
                    reducedS_[row][j] = normalize(reducedS_[row][j] - factor * reducedS_[current_row][j]);
                    rowOperation_[row][j] = normalize(rowOperation_[row][j] - factor * rowOperation_[current_row][j]);
                }
            }
        }
        
        current_row++;
    }
    
    // Compute dimNull_
    dimNull_ = n - rank;
    
    // Set z0_ (integer: 0 or 1)
    z0_ = (dimNull_ == 0) ? 1 : 0;

    // 计算分母和分子形式的rowOperation
    denoRowOperation_.resize(n);
    numeRowOperation_.resize(n, std::vector<PT>(n, PT(0)));
    
    for (int i = 0; i < n; ++i) {
        // 计算第i行所有元素分母的lcm
        denoRowOperation_[i] = lcmOfDenominators(rowOperation_[i]);
        
        // 计算分子形式: numeRowOperation_[i][j] = rowOperation_[i][j] * denoRowOperation_[i]
        for (int j = 0; j < n; ++j) {
            numeRowOperation_[i][j] = normalize(PT(rowOperation_[i][j] * RT(denoRowOperation_[i])));
        }
    }
}

// Find null space basis (for dimNull_ >= 1)
template<typename RT, typename PT>
std::vector<std::vector<RT>> Sector<RT, PT>::findNullSpace() const {
    if (dimNull_ == 0) {
        return std::vector<std::vector<RT>>(); // empty basis
    }
    
    int n = S_.size();
    std::vector<std::vector<RT>> null_basis;
    
    // Find free variables (columns without pivots in reducedS)
    std::vector<bool> is_pivot_col(n, false);
    std::vector<int> free_vars;
    
    // Mark pivot columns
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            if (!isZero(reducedS_[row][col])) {
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
        std::vector<RT> null_vector(n, RT(0));
        null_vector[free_var] = RT(1);
        
        // Back-substitute to find dependent variables
        for (int row = n - 1; row >= 0; --row) {
            // Find pivot column for this row
            int pivot_col = -1;
            for (int col = 0; col < n; ++col) {
                if (!isZero(reducedS_[row][col])) {
                    pivot_col = col;
                    break;
                }
            }
            
            if (pivot_col != -1 && pivot_col != free_var) {
                RT sum = RT(0);
                for (int col = pivot_col + 1; col < n; ++col) {
                    sum = sum + reducedS_[row][col] * null_vector[col];
                }
                null_vector[pivot_col] = RT(0) - sum;
            }
        }
        
        null_basis.push_back(null_vector);
    }
    
    return null_basis;
}

// Solve linear system (for dimNull_ == 0)
template<typename RT, typename PT>
std::vector<RT> Sector<RT, PT>::solveLinear(const std::vector<RT>& b) const {
    if (dimNull_ != 0) {
        throw std::runtime_error("Matrix is not invertible (dimNull_ != 0)");
    }
    
    // For dimNull_ == 0, rowOperation is the inverse matrix
    // Solve S * x = b => x = invS * b = rowOperation * b
    int n = S_.size();
    std::vector<RT> solution(n, RT(0));
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n && j < (int)b.size(); ++j) {
            solution[i] = solution[i] + rowOperation_[i][j] * b[j];
        }
        solution[i] = normalize(solution[i]);
    }
    
    return solution;
}

// Solve for c and z coefficients
template<typename RT, typename PT>
void Sector<RT, PT>::solveCandZ() {
    int n = S_.size(); // total matrix size (N+B)x(N+B)
    
    if (dimNull_ == 0) {
        // Case: dimNull_ == 0
        // S (C_1,...,C_B,z_1,...,z_N) = (1_1,...,1_B,0,...,0)
        std::vector<RT> rhs(n, RT(0));
        for (int i = 0; i < numBranch_ && i < n; ++i) {
            rhs[i] = RT(1);
        }
        candz_ = solveLinear(rhs);
    } else if (dimNull_ == 1) {
        // Case: dimNull_ == 1
        std::vector<std::vector<RT>> null_basis = findNullSpace();
        if (!null_basis.empty()) {
            candz_ = null_basis[0];
        }
    } else {
        // Case: dimNull_ >= 2
        std::vector<std::vector<RT>> null_basis = findNullSpace();
        RT best_C_sum = RT(0);
        std::vector<RT> best_solution(n, RT(0));
        bool found_nonzero = false;
        for (const auto& basis_vector : null_basis) {
            RT candidate_C_sum = RT(0);
            for (int i = 0; i < numBranch_ && i < (int)basis_vector.size(); ++i) {
                candidate_C_sum = candidate_C_sum + basis_vector[i];
            }
            if (!found_nonzero || !isZero(candidate_C_sum)) {
                best_C_sum = candidate_C_sum;
                best_solution = basis_vector;
                if (!isZero(candidate_C_sum)) {
                    found_nonzero = true;
                    break;
                }
            }
        }
        candz_ = best_solution;
    }
    
    // Normalize candz_
    for (int i = 0; i < (int)candz_.size(); ++i) {
        candz_[i] = normalize(candz_[i]);
    }
    
    // Compute C_ = sum of c_i
    C_ = RT(0);
    for (int i = 0; i < numBranch_; ++i) {
        C_ = C_ + candz_[i];
    }
    C_ = normalize(C_);

    // 计算分母和分子形式的candz
    denoCandZ_ = lcmOfDenominators(candz_);
    
    // 计算分子形式: numeCandZ_[i] = candz_[i] * denoCandZ_
    numeCandZ_.resize(candz_.size());
    for (size_t i = 0; i < candz_.size(); ++i) {
        numeCandZ_[i] = normalize(PT(candz_[i] * RT(denoCandZ_)));
    }
    
    // 计算 numeC_ = C_ * denoCandZ_
    numeC_ = normalize(PT(C_ * RT(denoCandZ_)));
}
