// 构造函数
template<typename T>
LinearSystem<T>::LinearSystem(const std::vector<std::vector<T>>& input_matrix) 
    : matrix(input_matrix), rows(input_matrix.size()), cols(input_matrix[0].size()) {}

// 执行高斯消元
template<typename T>
void LinearSystem<T>::eliminate() {
    TRACE("LinearSystem::eliminate");
    
    int current_row = 0;
    free_vars.clear();
    pivot_cols.resize(rows, -1);
    
    for (int col = 0; col < cols && current_row < rows; ++col) {
        // 找到当前列的非零元素
        int pivot_row = -1;
        for (int row = current_row; row < rows; ++row) {
            if (matrix[row][col] != T(0)) {
                pivot_row = row;
                break;
            }
        }
        
        // 如果当前列全为0，跳过
        if (pivot_row == -1) {
            continue;
        }
        
        // 交换行
        if (pivot_row != current_row) {
            std::swap(matrix[pivot_row], matrix[current_row]);
        }
        
        pivot_cols[current_row] = col;
        T pivot = matrix[current_row][col];
        
        // 将主元化为1
        for (int j = 0; j < cols; ++j) {
            matrix[current_row][j] = matrix[current_row][j] / pivot;
        }
        
        // 消元：将当前列的其他元素化为0
        for (int row = 0; row < rows; ++row) {
            if (row != current_row && matrix[row][col] != T(0)) {
                T factor = matrix[row][col];
                for (int j = 0; j < cols; ++j) {
                    matrix[row][j] = matrix[row][j] - factor * matrix[current_row][j];
                }
            }
        }
        
        current_row++;
    }
    
    // 找出自由变量
    std::vector<bool> is_pivot_col(cols, false);
    for (int i = 0; i < rows; ++i) {
        if (pivot_cols[i] != -1) {
            is_pivot_col[pivot_cols[i]] = true;
        }
    }
    
    for (int j = 0; j < cols; ++j) {
        if (!is_pivot_col[j]) {
            free_vars.push_back(j);
        }
    }
}

// 检查是否有非零解
template<typename T>
bool LinearSystem<T>::hasNontrivialSolution() {
    int rank = 0;
    for (int i = 0; i < rows; ++i) {
        if (pivot_cols[i] != -1) rank++;
    }
    return rank < cols; // 如果系数矩阵的秩小于变量个数，则有非零解
}

// 获取解空间的维度
template<typename T>
int LinearSystem<T>::getSolutionSpaceDimension() {
    int rank = 0;
    for (int i = 0; i < rows; ++i) {
        if (pivot_cols[i] != -1) rank++;
    }
    return cols - rank;
}

// 获取变量数量
template<typename T>
int LinearSystem<T>::getVariableCount() const {
    return cols;
}

template<typename T>
const std::vector<std::vector<T>>& LinearSystem<T>::getRREFMatrix() const {
    return matrix;
}

template<typename T>
const std::vector<int>& LinearSystem<T>::getPivotColumns() const {
    return pivot_cols;
}

template<typename T>
const std::vector<int>& LinearSystem<T>::getFreeVariableColumns() const {
    return free_vars;
}

// 输出解到流
template<typename T>
void LinearSystem<T>::printSolution(std::ostream& out, const std::vector<std::string>& var_names) {
    out << "{" << std::endl;
    
    // 所有变量的表达式
    std::vector<std::string> expressions(cols);
    
    // 初始化所有变量为自己
    for (int j = 0; j < cols; ++j) {
        expressions[j] = var_names[j];
    }
    
    // 从下往上回代，表达非自由变量
    for (int i = rows - 1; i >= 0; --i) {
        if (pivot_cols[i] == -1) continue; // 零行
        
        int pivot_col = pivot_cols[i];
        std::string expr = "";
        
        bool has_terms = false;
        
        for (int j = pivot_col + 1; j < cols; ++j) {
            if (matrix[i][j] != T(0)) {
                has_terms = true;
                T coeff = -matrix[i][j]; // 移项后变号
                
                if (!expr.empty()) {
                    expr += "+";
                }
                
                if (coeff == T(1)) {
                    expr += var_names[j];
                } else {
                    std::ostringstream oss;
                    oss << coeff;
                    expr += oss.str() + "*" + var_names[j];
                }
            }
        }
        
        if (!has_terms) {
            expressions[pivot_col] = "0";
        } else {
            expressions[pivot_col] = expr;
        }
    }
    
    // 输出所有变量的表达式
    for (int j = 0; j < cols; ++j) {
        out << var_names[j] << "->" << expressions[j];
        if (j < cols - 1) {
            out << ",";
        }
        out << std::endl;
    }
    
    out << "}" << std::endl;
}

// 解析线性系统矩阵文件
template<typename T>
std::vector<std::vector<T>> parseLinearSystem(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::string content((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();
    
    std::vector<std::vector<T>> matrix;
    size_t pos = 0;
    
    // 跳过开头的空白字符
    while (pos < content.size() && std::isspace(content[pos])) pos++;
    
    // 期望以 '{' 开始
    if (pos >= content.size() || content[pos] != '{') {
        throw std::runtime_error("Matrix should start with '{'");
    }
    pos++; // 跳过 '{'
    
    while (pos < content.size()) {
        // 跳过空白字符
        while (pos < content.size() && std::isspace(content[pos])) pos++;
        
        if (pos >= content.size()) break;
        
        if (content[pos] == '}') break; // 矩阵结束
        
        if (content[pos] == '{') {
            // 开始解析一行
            pos++; // 跳过 '{'
            std::vector<T> row;
            
            while (pos < content.size()) {
                // 跳过空白字符
                while (pos < content.size() && std::isspace(content[pos])) pos++;
                
                if (pos >= content.size()) break;
                if (content[pos] == '}') {
                    pos++; // 跳过 '}'
                    break;
                }
                
                // 解析数字
                if (content[pos] == ',' || content[pos] == '}') {
                    pos++;
                    continue;
                }
                
                std::string numStr;
                bool negative = false;
                if (content[pos] == '-') {
                    negative = true;
                    pos++;
                }
                
                while (pos < content.size() && std::isdigit(content[pos])) {
                    numStr += content[pos];
                    pos++;
                }
                
                if (!numStr.empty()) {
                    long long val = std::stoll(numStr);
                    if (negative) val = -val;
                    row.push_back(T(val));
                }
                
                // 跳过逗号
                while (pos < content.size() && (content[pos] == ',' || std::isspace(content[pos]))) pos++;
            }
            
            if (!row.empty()) {
                matrix.push_back(row);
            }
        } else {
            pos++; // 跳过其他字符
        }
        
        // 跳过行间的逗号
        while (pos < content.size() && (content[pos] == ',' || std::isspace(content[pos]))) pos++;
    }
    
    return matrix;
}

// 加载线性系统
template<typename T>
LinearSystem<T> loadLinearSystem(const std::string& filename) {
    TRACE("loadLinearSystem");
    
    // 解析矩阵数据
    auto matrix = parseLinearSystem<T>(filename);
    
    // 创建并返回线性系统对象
    LinearSystem<T> system(matrix);
    
    // 可以在这里添加其他初始化步骤，比如：
    // - 读取变量名（如果文件中包含）
    // - 验证矩阵格式
    // - 打印加载信息等
    
    return system;
}

// 读取变量名文件
std::vector<std::string> loadVariableNames(const std::string& filename) {
    TRACE("loadVariableNames");
    
    std::ifstream varsFile(filename);
    if (!varsFile) {
        throw std::runtime_error("Cannot open linearvars file: " + filename);
    }
    
    std::vector<std::string> var_names;
    std::string line;
    while (std::getline(varsFile, line)) {
        // 简单解析变量名，假设每行一个变量名
        if (!line.empty()) {
            var_names.push_back(line);
        }
    }
    varsFile.close();
    
    return var_names;
}
