#pragma once

#include "parser.hpp"
#include "diffeq.hpp"
#include "series.hpp"
#include "ffType.hpp"
#include "trace.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <atomic>
#include <thread>
#include <omp.h>
#include <stdexcept>
#include <sstream>

/// 积分级数计算配置
struct SeriesConfig {
    std::string work_dir;  // 工作目录，所有操作都在此目录下
    unsigned long long prime, a, b;
    int degree, num_threads;
};

/// 预计算的幂次缓存
template<typename T>
struct PowerCache {
    std::vector<T> neg_a_pows;
    std::vector<T> neg_b_pows;
    std::vector<T> one_minus_a_pows;
    std::vector<T> one_minus_b_pows;
    
    PowerCache(T a, T b, int degree);
};

/// 解析命令行参数并创建配置
SeriesConfig parseCommandLine(int argc, char* argv[]);

/// 设置计算环境（OpenMP、prime等）
template<typename T>
void setupEnvironment(const SeriesConfig& config);

/// 加载矩阵文件
template<typename T>
Matrix<T> loadMatrix(const std::string& filepath, Parser<T>& parser);

/// 验证矩阵维度的一致性
template<typename T>
void validateMatrices(const Matrix<T>& AX, 
                     const Matrix<T>& AY,
                     const Matrix<T>& red_coe);

/// 求解微分方程系统，返回主级数解 mseries[mi_bc][mi_red]
template<typename T>
std::vector<std::vector<Series<T>>> solveDiffEquations(
    const DiffSystem<T>& diffSys,
    int degree
);

/// 计算级数: series[intg][mi_bc] = sum_{mi_red} red_coe[intg][mi_red] * mseries[mi_bc][mi_red]
template<typename T>
std::vector<std::vector<Series<T>>> computeSeries(
    const Matrix<T>& red_coe,
    const std::vector<std::vector<Series<T>>>& mseries,
    int degree
);

/// 预计算积分权重
template<typename T>
PowerCache<T> precomputePowers(T a, T b, int degree);

/// 执行积分计算，将X^n Y^m级数转换为s的一元多项式
template<typename T>
std::vector<std::vector<std::vector<T>>> performIntegration(
    const std::vector<std::vector<Series<T>>>& series,
    const PowerCache<T>& powers,
    int degree
);

/// 输出积分级数结果
template<typename T>
void outputIntegratedSeries(
    const std::vector<std::vector<std::vector<T>>>& integrated_series,
    const SeriesConfig& config,
    int intg_count, int mi_count, int degree
);

/// 将单个积分级数格式化为字符串 (a[0]+a[1]s+a[2]s^2+...)
template<typename T>
std::string formatPolynomial(const std::vector<T>& coefficients, int degree);

#include "../src/integratedseries.tpp"
