#pragma once
#include <string>
#include <vector>

struct LoadConfig;
struct DHWConfig;

// Simple weather-driven load model: Q = base + UA*(Tin - Tout)
class LoadModel {
public:
    bool loadCSV(const std::string& path, const std::string& column);
    double demandAt(int hourIndex, const LoadConfig& cfg) const;
    double dhwAt(int hourIndex, const DHWConfig& cfg) const;
private:
    std::vector<double> toutC_; // outdoor temperature per step
};
