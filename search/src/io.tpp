static std::string searchTrim(const std::string& s) {
    size_t l = 0;
    while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) ++l;
    size_t r = s.size();
    while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) --r;
    return s.substr(l, r - l);
}

static std::vector<unsigned long long> parseSearchU64List(const std::string& raw) {
    std::string s = searchTrim(raw);
    size_t l = s.find('{');
    size_t r = s.rfind('}');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("List must be in braces: " + raw);
    }
    std::string body = s.substr(l + 1, r - l - 1);
    std::vector<unsigned long long> vals;
    std::stringstream ss(body);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = searchTrim(tok);
        if (tok.empty()) continue;
        vals.push_back(std::stoull(tok));
    }
    return vals;
}

SearchConfig parseSearchConfigFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::unordered_map<std::string, std::string> kv;
    std::string line;
    while (std::getline(in, line)) {
        line = searchTrim(line);
        if (line.empty()) continue;
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        kv[searchTrim(line.substr(0, eq))] = searchTrim(line.substr(eq + 1));
    }

    auto need = [&](const std::string& key) -> std::string {
        auto it = kv.find(key);
        if (it == kv.end()) {
            throw std::runtime_error("Missing key in config: " + key);
        }
        return it->second;
    };

    SearchConfig cfg;
    cfg.nuSize = std::stoi(need("N"));
    cfg.degreeD = std::stoi(need("deg"));
    cfg.numFBIMasters = static_cast<int>(parseSearchU64List(need("bc")).size());
    {
        auto it = kv.find("ncheck");
        if (it != kv.end()) {
            cfg.ncheck = std::stoi(it->second);
        }
    }
    if (cfg.ncheck < 0) {
        throw std::runtime_error("ncheck must be >= 0");
    }
    cfg.p = static_cast<mp_limb_t>(std::stoull(need("p")));
    return cfg;
}

int parseMaxSearchDegreeFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open maxsearchdeg file: " + path);
    }
    std::string line;
    while (std::getline(in, line)) {
        line = searchTrim(line);
        if (line.empty()) continue;
        return std::stoi(line);
    }
    throw std::runtime_error("maxsearchdeg file is empty: " + path);
}

std::vector<IntegralLabel> parseSearchTargetFile(const std::string& path, int expectedNuSize) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open target file: " + path);
    }

    std::vector<IntegralLabel> targets;
    std::string line;
    while (std::getline(in, line)) {
        line = searchTrim(line);
        if (line.empty()) continue;

        size_t l = line.find('{');
        size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid nu line in target file: " + line);
        }

        IntegralLabel label;
        std::string body = line.substr(l + 1, r - l - 1);
        std::stringstream ss(body);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = searchTrim(tok);
            if (!tok.empty()) {
                label.nu.push_back(std::stoi(tok));
            }
        }

        if (static_cast<int>(label.nu.size()) != expectedNuSize) {
            throw std::runtime_error("nu length mismatch in target file");
        }
        targets.push_back(std::move(label));
    }

    if (targets.empty()) {
        throw std::runtime_error("No nu entries in target file");
    }
    return targets;
}

template<typename T>
std::vector<SeriesSample<T>> parseSeriesFile(
    const std::string& path,
    const std::vector<IntegralLabel>& targets,
    int degreeD,
    int bcIndex) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open series file: " + path);
    }

    std::vector<SeriesSample<T>> samples;
    std::string line;
    int targetIndex = 0;
    while (std::getline(in, line)) {
        line = searchTrim(line);
        if (line.empty()) continue;
        if (targetIndex >= static_cast<int>(targets.size())) {
            throw std::runtime_error("Too many series lines in file: " + path);
        }

        size_t l = line.find('{');
        size_t r = line.rfind('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) {
            throw std::runtime_error("Invalid series line: " + line);
        }

        SeriesSample<T> sample;
        sample.label.integral = targets[targetIndex];
        sample.label.bcIndex = bcIndex;

        std::string body = line.substr(l + 1, r - l - 1);
        std::stringstream ss(body);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = searchTrim(tok);
            if (tok.empty()) continue;
            sample.coeffs.push_back(T(static_cast<unsigned long long>(std::stoull(tok))));
        }

        if (static_cast<int>(sample.coeffs.size()) != degreeD + 1) {
            throw std::runtime_error("Series degree mismatch in file: " + path);
        }

        samples.push_back(std::move(sample));
        ++targetIndex;
    }

    if (targetIndex != static_cast<int>(targets.size())) {
        throw std::runtime_error("Series line count mismatch in file: " + path);
    }
    return samples;
}
