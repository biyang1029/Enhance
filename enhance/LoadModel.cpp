#include "LoadModel.h"
#include "DataConfig.h"
#include <algorithm>
#include <fstream>
#include <sstream>

static inline std::string trim(const std::string& s){
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a==std::string::npos) return ""; return s.substr(a,b-a+1);
}

bool LoadModel::loadCSV(const std::string& path, const std::string& column){
    toutC_.clear();
    std::ifstream ifs(path.c_str());
    if (!ifs) return false;
    std::string line;
    // header
    if (!std::getline(ifs,line)) return false;
    std::vector<std::string> headers; std::stringstream ss(line); std::string cell;
    while(std::getline(ss,cell,',')) headers.push_back(trim(cell));
    int col = -1;
    for(size_t i=0;i<headers.size();++i){ if(headers[i]==column){ col=(int)i; break; } }
    if (col < 0){
        // assume first column is temperature and treat the header as first data row was read already
        std::stringstream ss2(line);
        if (std::getline(ss2,cell,',')){
            try { toutC_.push_back(std::stod(trim(cell))); } catch(...){}
        }
        col = 0;
    }
    while(std::getline(ifs,line)){
        if (line.empty()) continue;
        std::stringstream ls(line);
        // skip to target column
        for(int i=0;i<col; ++i){ if(!std::getline(ls,cell,',')) { cell.clear(); break; } }
        if (!std::getline(ls,cell,',')) continue;
        try { toutC_.push_back(std::stod(trim(cell))); } catch(...){}
    }
    return !toutC_.empty();
}

double LoadModel::demandAt(int hourIndex, const LoadConfig& cfg) const{
    if (cfg.enable_cutoff) {
        double ToutCut = cfg.indoor_T_C;
        if ( !toutC_.empty() && hourIndex >= 0){ ToutCut = toutC_[hourIndex % toutC_.size()]; } 
        if (ToutCut >= cfg.heat_cutoff_C) return 0.0;
    }
    if (!cfg.enable_weather) return cfg.base_kW;
    double Tout = cfg.indoor_T_C;
    if (!toutC_.empty() && hourIndex >= 0){
        Tout = toutC_[hourIndex % toutC_.size()];
    }
    double q = cfg.base_kW + cfg.UA_kW_per_K * (cfg.indoor_T_C - Tout);
    return std::max(0.0, q);
}

double LoadModel::dhwAt(int hourIndex, const DHWConfig& cfg) const{
    if (!cfg.enable) return 0.0;
    int h = 0;
    if (hourIndex >= 0) h = hourIndex % 24;
    double q = cfg.base_kW;
    auto in_window = [&](int start_h, int hours){
        if (hours <= 0) return false;
        int end = (start_h + hours - 1) % 24;
        if (hours >= 24) return true;
        if (start_h + hours <= 24) return (h >= start_h && h < start_h + hours);
        // wrap across midnight
        return (h >= start_h || h < end + 1);
    };
    if (in_window(cfg.morning_start_h, cfg.morning_hours)) q += std::max(0.0, cfg.morning_kW);
    if (in_window(cfg.evening_start_h, cfg.evening_hours)) q += std::max(0.0, cfg.evening_kW);
    return std::max(0.0, q);
}
