#include <iostream>
#include <stdexcept>
#include <string>

#include "../../../include/fi_pipeline.hpp"

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <S_path> <config_path> <target_path> <output_path>\n";
        return 1;
    }

    try {
        const std::string sPath = argv[1];
        const std::string configPath = argv[2];
        const std::string targetPath = argv[3];
        const std::string outputPath = argv[4];
        runFI1DSeriesPipeline(sPath, configPath, targetPath, outputPath);
        std::cout << "Pipeline finished.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
