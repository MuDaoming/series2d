#include "../include/integratedseries.hpp"
#include "../include/ffType.hpp"

int main(int argc, char* argv[]) {
    TRACE("main");
    
    try {
        // 解析命令行参数
        auto config = parseCommandLine(argc, argv);
        
        // 设置FlintMod模数
        FlintMod::set_modulus(config.prime);
        
        // 设置计算环境
        setupEnvironment<FlintMod>(config);
        
        using T = FlintMod;
        Parser<T> parser;
        
        // 加载矩阵文件
        auto AX = loadMatrix<T>(config.work_dir + "/DEX", parser);
        auto AY = loadMatrix<T>(config.work_dir + "/DEY", parser);
        auto red_coe = loadMatrix<T>(config.work_dir + "/coe", parser);
        
        // 验证矩阵维度
        validateMatrices(AX, AY, red_coe);
        
        int mi_count = AX.size();
        int intg_count = red_coe.size();
        
        // 创建微分方程系统
        DiffSystem<T> diffSys(AX, AY);
        
        // 求解微分方程
        auto mseries = solveDiffEquations(diffSys, config.degree);
        
        // 计算级数
        auto series = computeSeries(red_coe, mseries, config.degree);
        
        // 预计算积分权重
        auto powers = precomputePowers(T(config.a), T(config.b), config.degree);
        
        // 执行积分计算
        auto integrated_series = performIntegration(series, powers, config.degree);
        
        // 输出结果
        outputIntegratedSeries(integrated_series, config, intg_count, mi_count, config.degree);
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
