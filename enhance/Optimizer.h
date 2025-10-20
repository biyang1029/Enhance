#pragma once
#include "DataConfig.h"
#include <string>

class Optimizer {
public:
    explicit Optimizer(const DataConfig& cfg) : baseCfg_(cfg) {}
    // Returns true on success. Writes per-run outputs under runsRoot.
    bool run(const std::string& runsRoot);

private:
    DataConfig baseCfg_;
};

