#include <stdexcept>
#include <chrono>

// ===== DiffSystem模板类实现 =====

template<typename T>
DiffSystem<T>::DiffSystem(const std::vector<std::vector<Rational<T>>>& AX,
                          const std::vector<std::vector<Rational<T>>>& AY)
    : AX_(AX), AY_(AY) {
    validateSystem();
}

template<typename T>
void DiffSystem<T>::solve(std::vector<Series<T>>& result, const std::vector<T>& f0, int deg) {
    int nmi = getSystemSize();
    
    // 验证初值条件的大小
    if (f0.size() != static_cast<size_t>(nmi)) {
        throw std::invalid_argument("Initial condition vector size must match system size");
    }
    
    // 确保结果向量具有正确的大小
    result.resize(nmi);
    for (int k = 0; k < nmi; ++k) {
        result[k] = Series<T>(deg);
    }
    
    // 在循环外构造好耦合项，避免重复构造
    Series<T> g1(deg), g2(deg);
    
    // 依次求解每个方程
    for (int i = 0; i < nmi; ++i) {



        // 获取当前方程的系数（使用引用避免拷贝）
        const Rational<T>& R1 = AX_[i][i];  // X方向系数
        const Rational<T>& R2 = AY_[i][i];  // Y方向系数
        
        // 计算耦合项（使用引用传递避免拷贝）
        getSub(g1, g2, result, i);




    
        
        // 求解当前方程
        solveStandardDE(result[i], R1, g1, R2, g2, f0[i]);
    }
}

template<typename T>
void DiffSystem<T>::getSub(Series<T>& g1, Series<T>& g2, const std::vector<Series<T>>& f, int i) {
    // 从f[0]获取度数（如果f为空则使用默认度数0）
    int deg = (f.empty() || f[0].size() == 0) ? 0 : f[0].getDeg();
    
    // 确保g1和g2具有正确的度数
    if (g1.getDeg() != deg) g1 = Series<T>(deg);
    if (g2.getDeg() != deg) g2 = Series<T>(deg);
    
    // 如果是第一个方程，没有耦合项（g1和g2保持为零）
    if (i == 0) {
        return;
    }
    
    // 在循环外构造临时变量
    Series<T> contribution(deg);
    
    // 计算X方向的耦合项
    for (int j = 0; j < i; ++j) {
        if (AX_[i][j].isEmpty()) continue;
        
        // 使用静态函数避免临时对象构造
        // print deg for num and den of AX_[i][j]
;
        Series<T>::mulRat(contribution, f[j], AX_[i][j]);

        g1 += contribution;
    }
    g1 *= AX_[i][i].denominator;
    
    // 计算Y方向的耦合项
    for (int j = 0; j < i; ++j) {
        if (AY_[i][j].isEmpty()) continue;
        
        // 使用静态函数避免临时对象构造
        Series<T>::mulRat(contribution, f[j], AY_[i][j]);
        g2 += contribution;
    }
    g2 *= AY_[i][i].denominator;
}

template<typename T>
void DiffSystem<T>::solveStandardDE(Series<T>& result, const Rational<T>& R1, Series<T>& g1,
                                   const Rational<T>& R2, Series<T>& g2,
                                   const T& f0) {


    auto start = std::chrono::high_resolution_clock::now();
    // 从g1获取度数
    int deg = g1.getDeg();
    
    // 验证g1和g2的度数一致性
    if (g1.getDeg() != g2.getDeg()) {
        throw std::invalid_argument("g1 and g2 must have same degree");
    }
    
    // 确保result具有正确的度数
    if (result.getDeg() != deg) {
        throw std::invalid_argument("Result series degree does not match g1/g2 degree");
    }
    
    // 根据修改后的标准方程定义：
    // df/dX = R1 * f + g1
    // df/dY = R2 * f + g2
    // 需要将分母乘到方程左边，即：
    // Q1 * df/dX = P1 * f + Q1 * g1
    // Q2 * df/dY = P2 * f + Q2 * g2
    
    // // 预处理：将分母乘到g1和g2上（原地操作）
    // g1 *= R1.denominator;
    // g2 *= R2.denominator;
    
    // 获取分子分母多项式
    const auto& P1 = R1.numerator;
    const auto& Q1 = R1.denominator;
    const auto& P2 = R2.numerator;
    const auto& Q2 = R2.denominator;
    
    // 检查常数项是否为零
    T Q1_00 = Q1.getCoeff(0, 0);
    T Q2_00 = Q2.getCoeff(0, 0);
    if (Q1_00 == T(0) || Q2_00 == T(0)) {
        throw std::invalid_argument("Denominator constant term cannot be zero");
    }
    
    // 设置初值（其他系数默认为零）
    result.setCoeff(0, 0, f0);  // 设置初值
    
    // 获取各多项式的度数
    int degQ1 = Q1.getDegree();
    int degP1 = P1.getDegree();
    int degQ2 = Q2.getDegree();
    int degP2 = P2.getDegree();



    
    // 按照总度数递增的顺序计算系数
    for (int totalDeg = 1; totalDeg <= deg; ++totalDeg) {
        for (int m = 0; m <= totalDeg; ++m) {
            int n = totalDeg - m;
            
            if (n != 0) {
                // 使用X方向的微分方程
                // f[n,m] = (1/n) * (1/Q1[0,0]) * (g1[n-1,m] + Sum[P1[n-1-i,m-j]*f[i,j]] - Sum[(i+1)*Q1[n-1-i,m-j]*f[i+1,j]])
                
                T sum = g1.getCoeff(n-1, m);
                
                // 添加P1项
                int i_min = (n - 1 - degP1 > 0) ? (n - 1 - degP1) : 0;
                int i_max = (n - 1 < deg) ? (n - 1) : deg;
                int j_min = (m - degP1 > 0) ? (m - degP1) : 0;
                int j_max = (m < deg) ? m : deg;
                
                int numOfTerms = (i_max - i_min + 1) * (j_max - j_min + 1);
                int numOfMono = P1.getNumOfMonomials();
                // print numOfTerms and numOfMono

                if (numOfTerms < numOfMono)
                {   
                    // 循环的范围为 i,j
                    for (int i = i_min; i <= i_max; ++i) {
                        for (int j = j_min; j <= j_max; ++j) {
                            T P1_coeff = P1.getCoeff(n-1-i, m-j);
                            T f_coeff = result.getCoeff(i, j);
                            sum += P1_coeff * f_coeff;
                        }
                    }
                }
                else
                {
                    // 循环的范围为多项式的单项式
                    for (const auto& [power, coeff] : P1) {
                        int i = n - 1 - power.x_power;
                        int j = m - power.y_power;
                        
                        if (i >= 0 && j >= 0) {
                            T f_coeff = result.getCoeff(i, j);
                            sum += coeff * f_coeff;
                        }
                    }
                }
                
                // 减去Q1项
                i_min = (n - 1 - degQ1 > 0) ? (n - 1 - degQ1) : 0;
                i_max = (n - 1 < deg) ? (n - 1) : deg;
                j_min = (m - degQ1 > 0) ? (m - degQ1) : 0;
                j_max = (m < deg) ? m : deg;
                

                numOfTerms = (i_max - i_min + 1) * (j_max - j_min + 1);
                numOfMono = Q1.getNumOfMonomials();


                if (numOfTerms < numOfMono)
                {   
                    // 循环的范围为 i,j
                    for (int i = i_min; i <= i_max; ++i) {
                        for (int j = j_min; j <= j_max; ++j) {
                            if (i == n - 1 && j == m) continue; // 跳过( n-1, m )项
                            T Q1_coeff = Q1.getCoeff(n-1-i, m-j);
                            // T Q1_coeff = result.getCoeff(n-1-i, m-j);
                            T f_coeff = result.getCoeff(i+1, j);
                            sum -= T(i+1) * Q1_coeff * f_coeff;
                        }
                    }
                }
                else
                {
                    // 循环的范围为多项式的单项式
                    for (const auto& [power, coeff] : Q1) {
                        int i = n - 1 - power.x_power;
                        int j = m - power.y_power;
                        
                        if (i >= 0 && j >= 0) {
                            if (i == n - 1 && j == m) continue; // 跳过( n-1, m )项
                            T f_coeff = result.getCoeff(i+1, j);
                            sum -= T(i+1) * coeff * f_coeff;
                        }
                    }
                }
                result.setCoeff(n, m, sum / (T(n) * Q1_00));
                
            } else {
                // n == 0, 使用Y方向的微分方程
                // f[0,m] = (1/m) * (1/Q2[0,0]) * (g2[0,m-1] + Sum[P2[0-i,m-1-j]*f[i,j]] - Sum[(j+1)*Q2[0-i,m-1-j]*f[i,j+1]])
                
                T sum = g2.getCoeff(n, m-1);
                
                // 添加P2项
                int i_min = (n - degP2 > 0) ? (n - degP2) : 0;
                int i_max = (n < deg) ? n : deg;
                int j_min = (m - 1 - degP2 > 0) ? (m - 1 - degP2) : 0;
                int j_max = (m - 1 < deg) ? (m - 1) : deg;

                int numOfTerms = (i_max - i_min + 1) * (j_max - j_min + 1);
                int numOfMono = P2.getNumOfMonomials();

                // print numOfTerms and numOfMono

                if (numOfTerms < numOfMono)
                {
                    // 循环的范围为 i,j
                    for (int i = i_min; i <= i_max; ++i) {
                        for (int j = j_min; j <= j_max; ++j) {
                            T P2_coeff = P2.getCoeff(n-i, m-1-j);
                            T f_coeff = result.getCoeff(i, j);
                            sum += P2_coeff * f_coeff;
                        }
                    }
                }
                else
                {
                    // 循环的范围为多项式的单项式
                    for (const auto& [power, coeff] : P2) {
                        int i = n - power.x_power;
                        int j = m - 1 - power.y_power;

                        if (i >= 0 && j >= 0) {
                            T f_coeff = result.getCoeff(i, j);
                            sum += coeff * f_coeff;
                        }
                    }
                }
                
                // 减去Q2项
                i_min = (n - degQ2 > 0) ? (n - degQ2) : 0;
                i_max = (n < deg) ? n : deg;
                j_min = (m - 1 - degQ2 > 0) ? (m - 1 - degQ2) : 0;
                j_max = (m - 1 < deg) ? (m - 1) : deg;

                numOfTerms = (i_max - i_min + 1) * (j_max - j_min + 1);
                numOfMono = Q2.getNumOfMonomials();

                if (numOfTerms < numOfMono)
                {
                    // 循环的范围为 i,j
                    for (int i = i_min; i <= i_max; ++i) {
                        for (int j = j_min; j <= j_max; ++j) {
                            // 跳过(0, m-1)项
                            if (i == 0 && j == m-1) continue;
                            
                            T Q2_coeff = Q2.getCoeff(n-i, m-1-j);
                            T f_coeff = result.getCoeff(i, j+1);
                            sum -= T(j+1) * Q2_coeff * f_coeff;
                        }
                    }
                }
                else
                {
                    // 循环的范围为多项式的单项式
                    for (const auto& [power, coeff] : Q2) {
                        int i = n - power.x_power;
                        int j = m - 1 - power.y_power;
                        
                        // 跳过(0, m-1)项
                        if (i == n && j == m-1) continue;
                        
                        if (i >= 0 && j >= 0) {
                            T f_coeff = result.getCoeff(i, j+1);
                            sum -= T(j+1) * coeff * f_coeff;
                        }
                    }
                }
            
                
                result.setCoeff(0, m, sum / (T(m) * Q2_00));
            }
        }
    }


}


template<typename T>
int DiffSystem<T>::getSystemSize() const {
    return static_cast<int>(AX_.size());
}

template<typename T>
void DiffSystem<T>::validateSystem() const {
    if (AX_.size() != AY_.size()) {
        throw std::invalid_argument("AX and AY must have same size");
    }
    
    int n = static_cast<int>(AX_.size());
    for (int i = 0; i < n; ++i) {
        if (AX_[i].size() != static_cast<size_t>(n) || AY_[i].size() != static_cast<size_t>(n)) {
            throw std::invalid_argument("AX and AY must be square matrices");
        }
    }
}
