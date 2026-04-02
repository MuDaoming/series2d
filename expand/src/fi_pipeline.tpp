void runFI1DSeriesPipeline(const std::string& sPath,
                           const std::string& configPath,
                           const std::string& targetPath,
                           const std::string& outputPath) {
    namespace fs = std::filesystem;

    fs::path sP = fs::path(sPath);
    fs::path cfgP = fs::path(configPath);
    fs::path tgtP = fs::path(targetPath);
    fs::path outP = fs::path(outputPath);
    if (!sP.is_absolute()) sP = fs::current_path() / sP;
    if (!cfgP.is_absolute()) cfgP = fs::current_path() / cfgP;
    if (!tgtP.is_absolute()) tgtP = fs::current_path() / tgtP;
    if (!outP.is_absolute()) outP = fs::current_path() / outP;
    fs::path outFile = outP;
    if (!outFile.is_absolute()) outFile = fs::current_path() / outFile;
    if (outFile.has_parent_path()) {
        fs::create_directories(outFile.parent_path());
    }

    InputConfig cfg = parseConfigFile(cfgP.string());
    TargetConfig targetCfg = parseTargetFile(tgtP.string(), cfg.N);
    FlintMod::set_modulus(cfg.p);

    GiNaC::symbol X("X"), Y("Y");
    std::vector<std::vector<GiNaC::ex>> topS = parseMatrixFile(sP.string(), X, Y);
    if (topS.empty() || topS.size() != topS[0].size()) {
        throw std::runtime_error("S matrix must be non-empty and square");
    }
    if (static_cast<int>(topS.size()) != cfg.B + cfg.N) {
        throw std::runtime_error("S dimension mismatch: expected B+N");
    }

    const int numBranch = cfg.B;
    const int numProps = static_cast<int>(topS.size()) - numBranch;
    Family<GiNaC::ex, GiNaC::ex, GiNaC::ex> ginacFamily(topS, numProps, numBranch, X, Y);
    auto family = convertFamily(ginacFamily, X, Y);

    const FlintMod feynmanD(cfg.d);
    const FlintMod shiftA(cfg.a);
    const FlintMod shiftB(cfg.b);
    family.setMasterDelta(feynmanD);

    const int targetDeg = cfg.deg;
    SeriesSolver<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod> solver(family, targetDeg);
    if (cfg.bc.size() != static_cast<size_t>(solver.getNumMaster())) {
        std::ostringstream oss;
        oss << "config.bc size mismatch: got " << cfg.bc.size()
            << ", expected " << solver.getNumMaster() << ". ";
        oss << "master sectors: ";
        const auto& masterIdxs = family.getMasterIdxs();
        for (size_t i = 0; i < masterIdxs.size(); ++i) {
            std::vector<int> secvec = family.secvecFromIdx(masterIdxs[i]);
            oss << "{";
            for (size_t j = 0; j < secvec.size(); ++j) {
                oss << secvec[j];
                if (j + 1 < secvec.size()) oss << ",";
            }
            oss << "}";
            if (i + 1 < masterIdxs.size()) oss << ", ";
        }
        throw std::runtime_error(oss.str());
    }
    for (int i = 0; i < solver.getNumMaster(); ++i) {
        solver.setMasterBoundary(i, FlintMod(cfg.bc[i]));
    }
    solver.solve();

    IntegrandExpander<Rational<FlintMod>, Polynomial<FlintMod>, FlintMod> expander(
        solver, 2, targetDeg, feynmanD, shiftA, shiftB);
    IntegrationConfig<FlintMod> intCfg{shiftA, shiftB, targetDeg};
    SeriesIntegrator<FlintMod> integrator(intCfg);

    std::ofstream out(outFile.string());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output file: " + outFile.string());
    }

    for (const auto& nu : targetCfg.nus) {
        Series<FlintMod> fi2d = expander.getFI2DSeries(nu);
        std::vector<FlintMod> oneD = integrator.integrate(fi2d);

        out << "{";
        for (int d = 0; d <= targetDeg; ++d) {
            out << oneD[d].get_value();
            if (d < targetDeg) out << ",";
        }
        out << "}\n";
    }
}
