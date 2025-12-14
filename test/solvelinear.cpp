#include "../include/linear.hpp"
#include "../include/ffType.hpp"
#include "../include/trace.h"

int main(int argc, char* argv[]) {
    TRACE("main");
    
    try {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " <work_dir> <prime>" << std::endl;
            std::cerr << "Example: " << argv[0] << " . 998244353" << std::endl;
            return 1;
        }
        
        std::string work_dir = argv[1];
        unsigned long long prime = std::stoull(argv[2]);
        
        // 构建文件路径
        std::string matrix_file = work_dir + "/linearsys";
        std::string vars_file = work_dir + "/linearvars";
        std::string output_file = work_dir + "/linearsol";
        
        // 设置模数
        FlintMod::set_modulus(prime);
        
        using T = FlintMod;
        
        // 加载线性系统
        auto system = loadLinearSystem<T>(matrix_file);
        
        // 加载变量名
        auto var_names = loadVariableNames(vars_file);
        
        // 检查变量数量是否匹配
        if (var_names.size() != system.getVariableCount()) {
            throw std::runtime_error("Number of variables (" + std::to_string(var_names.size()) + 
                                   ") doesn't match matrix columns (" + std::to_string(system.getVariableCount()) + ")");
        }
        
        // 执行高斯消元
        system.eliminate();
        
        // 检查解的情况
        bool has_nontrivial = system.hasNontrivialSolution();
        
        if (has_nontrivial) {
            // 将解输出到文件
            TRACE_SCOPE("output_solution");
            std::ofstream output(output_file);
            if (!output) {
                throw std::runtime_error("Cannot create output file: " + output_file);
            }
            system.printSolution(output, var_names);
            output.close();
        } else {
            // 输出空解到文件
            TRACE_SCOPE("output_empty_solution");
            std::ofstream output(output_file);
            if (!output) {
                throw std::runtime_error("Cannot create output file: " + output_file);
            }
            output << "{}" << std::endl;
            output.close();
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
