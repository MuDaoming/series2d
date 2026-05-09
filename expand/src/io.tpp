static std::string trim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

static std::vector<mp_limb_t> parseU64List(const std::string& raw) {
    std::string s = trim(raw);
    size_t l = s.find('{');
    size_t r = s.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("List must be in braces: " + raw);
    }
    std::string body = s.substr(l + 1, r - l - 1);
    std::vector<mp_limb_t> vals;
    std::stringstream ss(body);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim(tok);
        if (tok.empty()) continue;
        vals.push_back(static_cast<mp_limb_t>(std::stoull(tok)));
    }
    return vals;
}

static std::vector<int> parseI32List(const std::string& raw) {
    std::string s = trim(raw);
    size_t l = s.find('{');
    size_t r = s.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("List must be in braces: " + raw);
    }
    std::string body = s.substr(l + 1, r - l - 1);
    std::vector<int> vals;
    std::stringstream ss(body);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim(tok);
        if (tok.empty()) continue;
        vals.push_back(std::stoi(tok));
    }
    return vals;
}

InputConfig parseConfigFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = trim(line.substr(0, eq));
        std::string val = trim(line.substr(eq + 1));
        kv[key] = val;
    }

    auto need = [&](const std::string& k) -> std::string {
        if (!kv.count(k)) throw std::runtime_error("Missing key in config: " + k);
        return kv[k];
    };

    InputConfig cfg;
    cfg.B = std::stoi(need("B"));
    cfg.N = std::stoi(need("N"));
    cfg.deg = std::stoi(need("deg"));
    cfg.p = static_cast<mp_limb_t>(std::stoull(need("p")));
    cfg.a = static_cast<mp_limb_t>(std::stoull(need("a")));
    cfg.b = static_cast<mp_limb_t>(std::stoull(need("b")));
    cfg.d = static_cast<mp_limb_t>(std::stoull(need("d")));
    if (kv.count("reduceMode")) {
        std::string mode = trim(kv["reduceMode"]);
        for (char& c : mode) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (mode == "maximalcut") mode = "maximal_cut";
        if (mode != "normal" && mode != "maximal_cut") {
            throw std::runtime_error("Invalid reduceMode in config: " + mode +
                                     " (expected normal or maximal_cut)");
        }
        cfg.reduceMode = mode;
    }
    if (kv.count("print2DMode")) {
        std::string mode = trim(kv["print2DMode"]);
        for (char& c : mode) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (mode != "target" && mode != "cache") {
            throw std::runtime_error("Invalid print2DMode in config: " + mode +
                                     " (expected target or cache)");
        }
        cfg.print2DMode = mode;
    }
    if (kv.count("sector")) {
        cfg.sector = parseI32List(kv["sector"]);
    }
    cfg.bc = parseU64List(need("bc"));
    return cfg;
}

TargetConfig parseTargetFile(const std::string& path, int expectedNuSize) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open target file: " + path);
    }

    TargetConfig targets;
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;

        // Backward compatible: ignore legacy "deg = ..." in target
        if (line.rfind("deg", 0) == 0) continue;

        size_t l = line.find('{');
        size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid nu line in target file: " + line);
        }

        std::string body = line.substr(l + 1, r - l - 1);
        std::vector<int> nu;
        std::stringstream ss(body);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = trim(tok);
            if (tok.empty()) continue;
            nu.push_back(std::stoi(tok));
        }
        if (static_cast<int>(nu.size()) != expectedNuSize) {
            throw std::runtime_error("nu length mismatch in target file");
        }
        targets.nus.push_back(std::move(nu));
    }

    if (targets.nus.empty()) {
        throw std::runtime_error("No nu entries in target file");
    }
    return targets;
}

// Keep this parser aligned with expand/test/SeriesSolver/test_series_solver.cpp
std::vector<std::vector<GiNaC::ex>> parseMatrixFile(const std::string& filename, GiNaC::symbol& X, GiNaC::symbol& Y) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    std::string processed;
    for (char c : content) {
        if (c != '\n' && c != '\r') {
            processed += c;
        }
    }

    size_t start = processed.find('{');
    size_t end = processed.rfind('}');
    if (start == std::string::npos || end == std::string::npos) {
        throw std::runtime_error("Invalid matrix format");
    }
    std::string matrixContent = processed.substr(start + 1, end - start - 1);

    std::vector<std::vector<GiNaC::ex>> matrix;
    int braceLevel = 0;
    std::string currentRow;
    for (size_t i = 0; i < matrixContent.size(); ++i) {
        char c = matrixContent[i];
        if (c == '{') {
            if (braceLevel == 0) {
                currentRow.clear();
            }
            braceLevel++;
            if (braceLevel > 1) currentRow += c;
        } else if (c == '}') {
            braceLevel--;
            if (braceLevel > 0) currentRow += c;
            if (braceLevel == 0 && !currentRow.empty()) {
                std::vector<GiNaC::ex> row;
                std::string element;
                int innerBrace = 0;
                for (char rc : currentRow) {
                    if (rc == '{') innerBrace++;
                    else if (rc == '}') innerBrace--;

                    if (rc == ',' && innerBrace == 0) {
                        GiNaC::symtab table;
                        table["X"] = X;
                        table["Y"] = Y;
                        GiNaC::parser reader(table);
                        try {
                            GiNaC::ex parsed = reader(element);
                            row.push_back(parsed);
                        } catch (...) {
                            row.push_back(GiNaC::ex(0));
                        }
                        element.clear();
                    } else {
                        element += rc;
                    }
                }
                if (!element.empty()) {
                    GiNaC::symtab table;
                    table["X"] = X;
                    table["Y"] = Y;
                    GiNaC::parser reader(table);
                    try {
                        GiNaC::ex parsed = reader(element);
                        row.push_back(parsed);
                    } catch (...) {
                        row.push_back(GiNaC::ex(0));
                    }
                }
                matrix.push_back(row);
            }
        } else if (braceLevel > 0) {
            currentRow += c;
        }
    }

    return matrix;
}
