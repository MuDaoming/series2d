#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <iterator>
#include "trace.h"

/// 模板化的线性系统类
template<typename T>
class LinearSystem {
private:
    std::vector<std::vector<T>> matrix;
    int rows, cols;
    std::vector<int> pivot_cols; // 记录每一行的主元列
    std::vector<int> free_vars;  // 自由变量的列索引
    
public:
    LinearSystem(const std::vector<std::vector<T>>& input_matrix);
    
    // 执行高斯消元
    void eliminate();
    
    // 检查是否有非零解
    bool hasNontrivialSolution();
    
    // 获取解空间的维度
    int getSolutionSpaceDimension();
    
    // 获取变量数量
    int getVariableCount() const;
    
    // 输出解到流
    void printSolution(std::ostream& out, const std::vector<std::string>& var_names);
};

/// 解析线性系统矩阵文件
template<typename T>
std::vector<std::vector<T>> parseLinearSystem(const std::string& filename);

/// 加载线性系统（包括矩阵和变量名）
template<typename T>
LinearSystem<T> loadLinearSystem(const std::string& filename);

/// 读取变量名文件
std::vector<std::string> loadVariableNames(const std::string& filename);

#include "../src/linear.tpp"
