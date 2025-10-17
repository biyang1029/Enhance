#pragma once
#include <string>
#include <vector>

struct LoadConfig;

// Simple weather-driven load model: Q = base + UA*(Tin - Tout)
class LoadModel {
public:
    bool loadCSV(const std::string& path, const std::string& column);
    double demandAt(int hourIndex, const LoadConfig& cfg) const;
private:
    std::vector<double> toutC_; // outdoor temperature per step
};

