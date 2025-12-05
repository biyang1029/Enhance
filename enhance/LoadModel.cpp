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
    // Detect EPW: skip first 8 header lines, use fixed dry-bulb column (index 6)
    auto ends_with = [](const std::string& s, const std::string& suf){
        if (s.size() < suf.size()) return false;
        return std::equal(suf.rbegin(), suf.rend(), s.rbegin(),
                          [](char a,char b){ return std::tolower(a)==std::tolower(b); });
    };
    bool isEpw = ends_with(path, ".epw");
    int col = 0;
    bool hasMDH = false; // year,month,day,hour,temp format (5+ columns)

    if (isEpw) {
        // Skip EPW header lines (8 lines)
        for (int i=0; i<8; ++i){ if (!std::getline(ifs,line)) return false; }
        col = 6; // dry-bulb temperature column in EPW (0-based index)
    } else {
        // CSV with header row
        if (!std::getline(ifs,line)) return false;
        std::vector<std::string> headers; std::stringstream ss(line); std::string cell;
        while(std::getline(ss,cell,',')) headers.push_back(trim(cell));
        col = -1;
        for(size_t i=0;i<headers.size();++i){ if(headers[i]==column){ col=(int)i; break; } }
        if (col < 0){
            // assume first column is temperature and treat the header as first data row was read already
            std::stringstream ss2(line);
            if (std::getline(ss2,cell,',')){
                try { toutC_.push_back(std::stod(trim(cell))); } catch(...){}
            }
            col = 0;
        }
    }

    // Parse loop
    struct Rec { int m=0,d=0,h=0; double t=0.0; };
    std::vector<Rec> recs;
    recs.reserve(9000);

    while(std::getline(ifs,line)){
        if (line.empty()) continue;
        std::stringstream ls(line);
        std::string cell;
        std::vector<std::string> cells;
        while(std::getline(ls,cell,',')) cells.push_back(trim(cell));
        if (cells.empty()) continue;
        if (cells.size() >= 5) {
            // year, month, day, hour(1-24), temperature (col 4 index)
            try {
                Rec r;
                r.m = std::stoi(cells[1]);
                r.d = std::stoi(cells[2]);
                r.h = std::stoi(cells[3]);
                r.t = std::stod(cells[4]);
                recs.push_back(r);
                hasMDH = true;
                continue;
            } catch(...) {
                // fall through to generic column handling
            }
        }
        // generic column extraction
        int targetCol = col;
        if ((int)cells.size() <= targetCol) continue;
        try { toutC_.push_back(std::stod(cells[targetCol])); } catch(...) {}
    }

    if (hasMDH && !recs.empty()) {
        // rotate to start at first 10/15 hour=1 if found
        size_t startIdx = 0;
        for (size_t i=0;i<recs.size();++i){
            if (recs[i].m==10 && recs[i].d==15){ startIdx = i; break; }
        }
        for (size_t k=0;k<recs.size();++k){
            const Rec& r = recs[(startIdx + k) % recs.size()];
            toutC_.push_back(r.t);
        }
    }

    return !toutC_.empty();
}

double LoadModel::demandAt(int hourIndex, const LoadConfig& cfg) const{
    if (cfg.enable_cutoff) {
        double ToutCut = toutAt(hourIndex, cfg.indoor_T_C);
        if (ToutCut >= cfg.heat_cutoff_C) return 0.0;
    }
    if (!cfg.enable_weather) return cfg.base_kW;
    double Tout = toutAt(hourIndex, cfg.indoor_T_C);
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

double LoadModel::toutAt(int hourIndex, double fallback) const {
    if (!toutC_.empty() && hourIndex >= 0) {
        return toutC_[hourIndex % toutC_.size()];
    }
    return fallback;
}
