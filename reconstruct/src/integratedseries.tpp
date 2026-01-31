// PowerCache 构造函数实现
template<typename T>
PowerCache<T>::PowerCache(T a, T b, int degree) {
    T minusa = -a;
    T minusb = -b;
    T oneminusa = T(1) - a;
    T oneminusb = T(1) - b;
    
    // 初始化幂次向量
    neg_a_pows.resize(degree + 2, T(1));
    neg_b_pows.resize(degree + 2, T(1));
    one_minus_a_pows.resize(degree + 2, T(1));
    one_minus_b_pows.resize(degree + 2, T(1));
    
    // 预计算所有幂次
    for (int k = 1; k <= degree + 1; ++k) {
        neg_a_pows[k] = neg_a_pows[k - 1] * minusa;
        neg_b_pows[k] = neg_b_pows[k - 1] * minusb;
        one_minus_a_pows[k] = one_minus_a_pows[k - 1] * oneminusa;
        one_minus_b_pows[k] = one_minus_b_pows[k - 1] * oneminusb;
    }
}

// 解析命令行参数并创建配置
SeriesConfig parseCommandLine(int argc, char* argv[]) {
    if (argc != 6) {
        throw std::runtime_error("Usage: <program> <work_dir> <p> <a> <b> <deg>");
    }
    
    SeriesConfig config;
    config.work_dir = argv[1];
    config.prime = std::stoull(argv[2]);
    config.a = std::stoull(argv[3]);
    config.b = std::stoull(argv[4]);
    config.degree = std::atoi(argv[5]);
    config.num_threads = 6;  // 默认线程数
    
    return config;
}

// 设置计算环境
template<typename T>
void setupEnvironment(const SeriesConfig& config) {
    // 设置OpenMP线程数
    // omp_set_num_threads(config.num_threads);
    
    // 设置质数（这里假设T是FlintMod）
    if constexpr (std::is_same_v<T, FlintMod>) {
        FlintMod::set_modulus(config.prime);
    }
}

// 加载矩阵文件
template<typename T>
Matrix<T> loadMatrix(const std::string& filepath, Parser<T>& parser) {
    TRACE("loadMatrix");
    
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Error: Cannot open file " + filepath);
    }
    
    std::string content((std::istreambuf_iterator<char>(file)),
                       std::istreambuf_iterator<char>());
    file.close();
    
    return parser.parseMatrix(content);
}

// 验证矩阵维度的一致性
template<typename T>
void validateMatrices(const Matrix<T>& AX, 
                     const Matrix<T>& AY,
                     const Matrix<T>& red_coe) {
    // 验证AX和AY矩阵维度
    if (AX.size() != AY.size() || AX[0].size() != AY[0].size()) {
        throw std::runtime_error("Error: AX and AY matrices must have same dimensions");
    }
    
    int mi_count = AX.size();
    
    // 验证约化系数矩阵
    int intg_count = red_coe.size();
    int j_count = red_coe[0].size();
    
    if (j_count != mi_count) {
        throw std::runtime_error("Error: Reduction coefficient matrix columns (" + 
                                std::to_string(j_count) + ") must equal system size (" + 
                                std::to_string(mi_count) + ")");
    }
}

// 求解微分方程系统
template<typename T>
std::vector<std::vector<Series<T>>> solveDiffEquations(
    const DiffSystem<T>& diffSys,
    int degree) {
    
    TRACE("solveDiffEquations");
    
    int mi_count = diffSys.getSystemSize();
    
    // 计算 mseries[mi_bc][mi_red] - 对每个边界条件求解
    // std::vector<std::vector<Series<T>>> mseries(mi_count, std::vector<Series<T>>(mi_count, Series<T>(degree)));

    
    // debug
    std::vector<std::vector<Series<T>>> mseries(1, std::vector<Series<T>>(mi_count, Series<T>(degree)));


        
    // 使用OpenMP并行化边界条件求解
    // #pragma omp parallel for 
    // for (int mi_bc = 0; mi_bc < mi_count; ++mi_bc) {
    // debug
    for (int mi_bc = 0; mi_bc < 1; ++mi_bc) {
        // mi_bc标记边界条件
        std::vector<T> boundary_conditions(mi_count, T(0));
        boundary_conditions[mi_bc] = T(1);  // 边界条件第mi_bc个分量为1
        
        // 求解微分方程
        // 注意：需要使用const_cast来调用非const的solve方法
        const_cast<DiffSystem<T>&>(diffSys).solve(mseries[mi_bc], boundary_conditions, degree);
    }
    
    return mseries;
}

// 计算级数
template<typename T>
std::vector<std::vector<Series<T>>> computeSeries(
    const Matrix<T>& red_coe,
    const std::vector<std::vector<Series<T>>>& mseries,
    int degree) {
    
    TRACE("computeSeries");
    
    int intg_count = red_coe.size();
    // int mi_count = mseries.size();
    int mi_count = mseries[0].size();
    
    // 计算 series[intg][mi_bc] = sum_{mi_red} red_coe[intg][mi_red] * mseries[mi_bc][mi_red]
    // std::vector<std::vector<Series<T>>> series(intg_count, std::vector<Series<T>>(mi_count, Series<T>(degree)));

    // debug
    std::vector<std::vector<Series<T>>> series(intg_count, std::vector<Series<T>>(1, Series<T>(degree)));

    // 进度监控, 粒度为 (intg, mi_bc), 对齐OpenMP并行化
    // std::atomic<size_t> completed{0};
    // // size_t total_tasks = static_cast<size_t>(intg_count) * static_cast<size_t>(mi_count);
    // // debug
    // size_t total_tasks = static_cast<size_t>(intg_count) * static_cast<size_t>(1);
    // std::atomic<bool> stopMonitor{false};

    // 监控线程：每5秒打印已完成任务 / 总任务
    // std::thread monitor([&]() {
    //     while (!stopMonitor.load(std::memory_order_relaxed)) {
    //         size_t cur = completed.load(std::memory_order_relaxed);
    //         std::cout << "Progress: " << cur << "/" << total_tasks << std::endl;
    //         if (cur >= total_tasks) break;
    //         std::this_thread::sleep_for(std::chrono::seconds(5));
    //     }
    // });
    // 使用OpenMP并行化外层两个循环（collapse)
    // #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int intg = 0; intg < intg_count; ++intg) {
        // for (int mi_bc = 0; mi_bc < mi_count; ++mi_bc) {
        // debug
        for (int mi_bc = 0; mi_bc < 1; ++mi_bc) {
            // series[intg][mi_bc] = sum_{mi_red} red_coe[intg][mi_red] * mseries[mi_bc][mi_red]
            Series<T> temp(degree);
            temp.setZero();
            for (int mi_red = 0; mi_red < mi_count; ++mi_red) {
                // temp += red_coe[intg][mi_red] * mseries[mi_bc][mi_red]
                Series<T> product(degree);
                Series<T>::mulRat(product, mseries[mi_bc][mi_red], red_coe[intg][mi_red]);
                temp += product;
            }
            series[intg][mi_bc] = temp;
            // 每完成一个 (intg, mi_bc) 迭代即算一个任务完成
            // completed.fetch_add(1, std::memory_order_relaxed);
        }
    }

    // 通知监控线程结束并 join
    // stopMonitor.store(true, std::memory_order_relaxed);
    // if (monitor.joinable()) monitor.join();

    // 最终确保打印 100%
    // std::cout << "Progress: " << total_tasks << "/" << total_tasks << std::endl;

    return series;
}

// 预计算积分权重
template<typename T>
PowerCache<T> precomputePowers(T a, T b, int degree) {
    return PowerCache<T>(a, b, degree);
}

// 执行积分计算
template<typename T>
std::vector<std::vector<std::vector<T>>> performIntegration(
    const std::vector<std::vector<Series<T>>>& series,
    const PowerCache<T>& powers,
    int degree) {
    
    TRACE("performIntegration");
    
    int intg_count = series.size();
    int mi_count = series[0].size();
    
    
    // integrated_series 有维度 intg_count * mi_count * (degree+1)
    // std::vector<std::vector<std::vector<T>>> integrated_series(
    //     intg_count, std::vector<std::vector<T>>(mi_count, std::vector<T>(degree + 1)));
    // debug
    std::vector<std::vector<std::vector<T>>> integrated_series(
    intg_count, std::vector<std::vector<T>>(1, std::vector<T>(degree + 1)));
    
    // 积分权重计算 - 使用OpenMP并行化前三层循环
    // #pragma omp parallel for collapse(3)
    for (int intg = 0; intg < intg_count; ++intg) {
        // for (int mi_bc = 0; mi_bc < mi_count; ++mi_bc) {
        // debug
        for (int mi_bc = 0; mi_bc < 1; ++mi_bc) {
            for (int d = 0; d <= degree; ++d) {
                integrated_series[intg][mi_bc][d] = T(0);
                for (int n_idx = 0; n_idx <= d; ++n_idx) {
                    int m_idx = d - n_idx;
                    T coeff = series[intg][mi_bc].getCoeff(n_idx, m_idx);
                    if (coeff == T(0)) continue;

                    // compute integration weight
                    T weight_n = (powers.one_minus_a_pows[n_idx + 1] - powers.neg_a_pows[n_idx + 1]) / T(n_idx + 1);
                    T weight_m = (powers.one_minus_b_pows[m_idx + 1] - powers.neg_b_pows[m_idx + 1]) / T(m_idx + 1);
                    T total_weight = weight_n * weight_m;

                    // apply weight
                    integrated_series[intg][mi_bc][d] += coeff * total_weight;
                }
            }
        }
    }
    
    return integrated_series;
}

// 输出积分级数结果
template<typename T>
void outputIntegratedSeries(
    const std::vector<std::vector<std::vector<T>>>& integrated_series,
    const SeriesConfig& config,
    int intg_count, int mi_count, int degree) {
    
    // output each integrated_series[intg][mi_bc] as a[0]+a[1]s+a[2]s^2+...+a[degree]s^degree
    std::ofstream integratedOutFile(config.work_dir + "/series");
    if (!integratedOutFile) {
        throw std::runtime_error("Error: Cannot create integrated output file");
    }
    
    integratedOutFile << "{";
    for (int intg = 0; intg < intg_count; ++intg) {
        if (intg > 0) integratedOutFile << ",";
        integratedOutFile << "{";
        // for (int mi_bc = 0; mi_bc < mi_count; ++mi_bc) {
        // debug
        for (int mi_bc = 0; mi_bc < 1; ++mi_bc) {
            if (mi_bc > 0) integratedOutFile << ",";
            
            // output integrated_series[intg][mi_bc] as a[0]+a[1]s+a[2]s^2+...+a[degree]s^degree
            integratedOutFile << formatPolynomial(integrated_series[intg][mi_bc], degree);
        }
        integratedOutFile << "}";
    }
    integratedOutFile << "}" << std::endl;
    integratedOutFile.close();
}

// 将单个积分级数格式化为字符串
template<typename T>
std::string formatPolynomial(const std::vector<T>& coefficients, int degree) {
    std::ostringstream oss;
    bool first_term = true;
    
    for (int d = 0; d <= degree; ++d) {
        T coeff = coefficients[d];
        if (coeff == T(0)) continue; // 跳过零系数项
        
        if (!first_term) {
            oss << "+";
        }
        first_term = false;
        
        if (d == 0) {
            // 常数项
            oss << coeff;
        } else if (d == 1) {
            // 一次项
            if (coeff == T(1)) {
                oss << "s";
            } else {
                oss << coeff << "*s";
            }
        } else {
            // 高次项
            if (coeff == T(1)) {
                oss << "s^" << d;
            } else {
                oss << coeff << "*s^" << d;
            }
        }
    }
    
    // 如果所有系数都为0，输出0
    if (first_term) {
        oss << "0";
    }
    
    return oss.str();
}
