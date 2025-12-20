// Constructor
template<typename T>
Sector<T>::Sector(const std::vector<std::vector<T>>& matrix, int n_props, int n_branch) 
    : S(matrix), numProps(n_props), numBranch(n_branch), dimNull(-1), z0(T(0)), C(T(0)) {
    
    // Initialize candz vector
    candz.resize(numBranch + numProps);
    
    // Perform row reduction
    rowReduce();
    
    // Solve for c and z
    solveCandZ();
}

// Perform RREF and compute rowOperation matrix simultaneously
template<typename T>
void Sector<T>::rowReduce() {
    int n = S.size();
    
    // Initialize reducedS as copy of S
    reducedS = S;
    
    // Initialize rowOperation as identity matrix
    rowOperation.assign(n, std::vector<T>(n, T(0)));
    for (int i = 0; i < n; ++i) {
        rowOperation[i][i] = T(1);
    }
    
    int rank = 0;
    int current_row = 0;
    
    // Perform RREF with simultaneous operations on rowOperation
    for (int col = 0; col < n && current_row < n; ++col) {
        // Find pivot in current column
        int pivot_row = -1;
        for (int row = current_row; row < n; ++row) {
            if (reducedS[row][col] != T(0)) {
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
            std::swap(reducedS[pivot_row], reducedS[current_row]);
            std::swap(rowOperation[pivot_row], rowOperation[current_row]);
        }
        
        T pivot = reducedS[current_row][col];
        
        // Scale pivot row to make pivot = 1
        for (int j = 0; j < n; ++j) {
            reducedS[current_row][j] = reducedS[current_row][j] / pivot;
            rowOperation[current_row][j] = rowOperation[current_row][j] / pivot;
        }
        
        // Eliminate other rows (both above and below)
        for (int row = 0; row < n; ++row) {
            if (row != current_row && reducedS[row][col] != T(0)) {
                T factor = reducedS[row][col];
                for (int j = 0; j < n; ++j) {
                    reducedS[row][j] = reducedS[row][j] - factor * reducedS[current_row][j];
                    rowOperation[row][j] = rowOperation[row][j] - factor * rowOperation[current_row][j];
                }
            }
        }
        
        current_row++;
    }
    
    // Compute dimNull
    // int rank = 0;
    // for (int row = 0; row < n; ++row) {
    //     bool is_zero_row = true;
    //     for (int col = 0; col < n; ++col) {
    //         if (reducedS[row][col] != T(0)) {
    //             is_zero_row = false;
    //             break;
    //         }
    //     }
    //     if (!is_zero_row) {
    //         rank++;
    //     }
    // }
    dimNull = n - rank;
    
    // Set z0
    z0 = (dimNull == 0) ? T(1) : T(0);
}

// Find null space basis (for dimNull >= 1)
template<typename T>
std::vector<std::vector<T>> Sector<T>::findNullSpace() const {
    if (dimNull == 0) {
        return std::vector<std::vector<T>>(); // empty basis
    }
    
    int n = S.size();
    std::vector<std::vector<T>> null_basis;
    
    // Find free variables (columns without pivots in reducedS)
    std::vector<bool> is_pivot_col(n, false);
    std::vector<int> free_vars;
    
    // Mark pivot columns
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            if (reducedS[row][col] != T(0)) {
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
                if (reducedS[row][col] != T(0)) {
                    pivot_col = col;
                    break;
                }
            }
            
            if (pivot_col != -1 && pivot_col != free_var) {
                T sum = T(0);
                for (int col = pivot_col + 1; col < n; ++col) {
                    sum += reducedS[row][col] * null_vector[col];
                }
                null_vector[pivot_col] = -sum;
            }
        }
        
        null_basis.push_back(null_vector);
    }
    
    return null_basis;
}

// Solve linear system (for dimNull == 0)
template<typename T>
std::vector<T> Sector<T>::solveLinear(const std::vector<T>& b) const {
    if (dimNull != 0) {
        throw std::runtime_error("Matrix is not invertible (dimNull != 0)");
    }
    
    // For dimNull == 0, rowOperation is the inverse matrix
    // Solve S * x = b => x = invS * b = rowOperation * b
    int n = S.size();
    std::vector<T> solution(n, T(0));
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n && j < b.size(); ++j) {
            solution[i] += rowOperation[i][j] * b[j];
        }
    }
    
    return solution;
}

// Solve for c and z coefficients
template<typename T>
void Sector<T>::solveCandZ() {
    int n = S.size(); // total matrix size (N+B)x(N+B)
    
    if (dimNull == 0) {
        // Case: dimNull == 0
        // S (C_1,...,C_B,z_1,...,z_N) = (1_1,...,1_B,0,...,0)
        std::vector<T> rhs(n, T(0));
        for (int i = 0; i < numBranch && i < n; ++i) {
            rhs[i] = T(1);
        }
        candz = solveLinear(rhs);
    } else if (dimNull == 1) {
        // Case: dimNull == 1
        std::vector<std::vector<T>> null_basis = findNullSpace();
        if (!null_basis.empty()) {
            candz = null_basis[0];
        }
    } else {
        // Case: dimNull >= 2
        std::vector<std::vector<T>> null_basis = findNullSpace();
        T best_C_sum = T(0);
        std::vector<T> best_solution(n, T(0));
        bool found_nonzero = false;
        for (const auto& basis_vector : null_basis) {
            T candidate_C_sum = T(0);
            for (int i = 0; i < numBranch && i < basis_vector.size(); ++i) {
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
        candz = best_solution;
    }
    
    // Compute C = sum of c_i
    C = T(0);
    for (int i = 0; i < numBranch; ++i) {
        C += candz[i];
    }
}

template<typename T>
std::vector<std::vector<int>> Sector<T>::getBranch() const {
    std::vector<std::vector<int>> branch(numBranch);

    for (int i = 0; i < numBranch; i++) {
        for (int j = numBranch; j < numBranch + numProps; j++) {
            if (S[i][j] == T(1))
                branch[i].push_back(j - numBranch);
        }
    }

    return branch;
}

// 判断删除指定的传播子后，是否有某个branch的所有元素都被删除
template<typename T>
bool Sector<T>::isDeleteBranch(const std::vector<std::vector<int>>& branch, const std::vector<int>& props) {
    // 遍历每个 branch
    for (const auto& br : branch) {
        // 如果这个 branch 为空，跳过
        if (br.empty()) {
            continue;
        }
        
        // 检查这个 branch 的所有元素是否都在 props 中
        bool all_deleted = true;
        for (int prop_idx : br) {
            // 检查 prop_idx 是否在 props 中
            bool found = false;
            for (int p : props) {
                if (p == prop_idx) {
                    found = true;
                    break;
                }
            }
            // 如果有一个元素不在 props 中，说明这个 branch 不会被完全删除
            if (!found) {
                all_deleted = false;
                break;
            }
        }
        
        // 如果这个 branch 的所有元素都在 props 中，返回 true
        if (all_deleted) {
            return true;
        }
    }
    
    // 没有任何 branch 被完全删除
    return false;
}

// 从 topS 和 nu 得到 subS
template<typename T>
std::vector<std::vector<T>> Sector<T>::getSubS(const std::vector<std::vector<T>>& topS, 
                                                 const std::vector<int>& nu, 
                                                 int numBranch) {
    int n = topS.size();
    int numProps = nu.size();
    
    // 计算 subS 的维度：B + nu 中为 1 的元素个数
    int sub_n = numBranch;
    for (int i = 0; i < numProps; i++) {
        if (nu[i] == 1) {
            sub_n++;
        }
    }
    
    // 创建 subS
    std::vector<std::vector<T>> subS(sub_n, std::vector<T>(sub_n));
    
    int subrow = 0;
    for (int i = 0; i < n; i++) {
        // 只选择前 B 行或者 nu[i-B] == 1 的行
        if (i >= numBranch && nu[i - numBranch] != 1) continue;
        
        int subcol = 0;
        for (int j = 0; j < n; j++) {
            // 只选择前 B 列或者 nu[j-B] == 1 的列
            if (j >= numBranch && nu[j - numBranch] != 1) continue;
            
            subS[subrow][subcol] = topS[i][j];
            subcol++;
        }
        subrow++;
    }
    
    return subS;
}
