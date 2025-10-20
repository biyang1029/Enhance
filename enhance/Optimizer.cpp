#define _CRT_SECURE_NO_WARNINGS 1
#include "Optimizer.h"
#include "SimulationController.h"
#include <algorithm>
#include <direct.h>
#include <fstream>
#include <sstream>

using std::string;
// simple directory creation helper (recursive, using _mkdir)
static void ensure_dir(const std::string& path){
    if (path.empty()) return;
    std::string norm; norm.reserve(path.size());
    for(char c : path) norm.push_back(c=='\\'?'/':c);
    std::string cur;
    for(size_t i=0;i<norm.size();++i){
        char c = norm[i]; cur.push_back(c);
        if (c=='/' || i==norm.size()-1){
            if (!cur.empty() && (cur.back()=='/')){
                if (cur.size()>1) _mkdir(cur.c_str());
            } else {
                _mkdir(cur.c_str());
            }
        }
    }
}

static inline const char* getenv_c(const char* k){ const char* v = std::getenv(k); return v ? v : ""; }

static std::vector<double> parseLengthGrid(const std::string& spec){
    std::vector<double> vals;
    if (spec.find(':') != std::string::npos){
        double a=0,b=0; int n=0; char c1=':', c2=':'; std::stringstream ss(spec); string s;
        // simple parse: a:b:n
        size_t p1 = spec.find(':'), p2 = spec.find(':', p1+1);
        try {
            a = std::stod(spec.substr(0,p1));
            b = std::stod(spec.substr(p1+1, p2-p1-1));
            n = std::stoi(spec.substr(p2+1));
        } catch(...) { n=0; }
        if (n <= 1){ vals.push_back(a); }
        else {
            double step = (b - a) / double(n-1);
            for (int i=0;i<n;++i) vals.push_back(a + i*step);
        }
    } else {
        std::stringstream ss(spec); string tok;
        while(std::getline(ss, tok, ',')){
            try { if(!tok.empty()) vals.push_back(std::stod(tok)); } catch(...) {}
        }
    }
    return vals;
}

static bool readSummary(const std::string& summaryCsv,
                        double& comp_kWh, double& pump_kWh, double& heat_kWh, double& LCOH_annual){
    std::ifstream ifs(summaryCsv.c_str()); if(!ifs) return false;
    string line; if (!std::getline(ifs, line)) return false; // header
    if (!std::getline(ifs, line)) return false;
    std::stringstream ss(line);
    // header is: elec_kWh,comp_kWh,pump_kWh,heat_kWh,space_kWh,dhw_kWh,cost_elec,capex,lifetime_y,discount,CRF,LCOH_annual
    // we need columns 2,3,4,12 (1-based)
    std::vector<string> cells; string cell; while(std::getline(ss, cell, ',')) cells.push_back(cell);
    auto rd = [&](int idx)->double{ try{ return std::stod(cells.at(idx)); } catch(...) { return 0.0; } };
    comp_kWh = rd(1); pump_kWh = rd(2); heat_kWh = rd(3); LCOH_annual = rd(11);
    return true;
}

static double computeUnmetFromResults(const std::string& resultsCsv, double dt_h){
    std::ifstream ifs(resultsCsv.c_str()); if(!ifs) return 0.0;
    string line; if (!std::getline(ifs, line)) return 0.0; // header
    double sum=0.0; while(std::getline(ifs, line)){
        if(line.empty()) continue; std::stringstream ss(line); std::vector<string> cells; string c; while(std::getline(ss,c,',')) cells.push_back(c);
        int idxUnmet = 14; // step,T_source,T_return,COP,Q_out,P_el,Q_geo,model,T_tank,HP_on,Q_space_req,Q_dhw_req,Q_space_serv,Q_dhw_serv,Q_unmet,... => index 14
        if ((int)cells.size() > idxUnmet){ try { sum += std::stod(cells[idxUnmet]) * dt_h; } catch(...){} }
    }
    return sum;
}

bool Optimizer::run(const std::string& runsRoot){
    ensure_dir(runsRoot);
    // Length grid from env LEN_GRID; default 0 .. 30% depth with 7 points
    double depth = baseCfg_.well.depth_m;
    string lenSpec = getenv_c("LEN_GRID");
    if (lenSpec.empty()){
        double Lmax = 0.3 * depth; std::ostringstream os; os<<0<<":"<<Lmax<<":"<<7; lenSpec = os.str();
    }
    auto Lvals = parseLengthGrid(lenSpec);
    if (Lvals.empty()) Lvals.push_back(0.0);

    // For each value, run simulation
    _putenv_s("SUMMARY_ONLY", "1");
    const double dt_h = baseCfg_.time.timeStep_s / 3600.0;
    struct Row{ double L_enh_m,z0,z1,comp,pump,heat,elec,scop,LCOH,unmet; };
    std::vector<Row> rows;
    for (double Lm : Lvals){
        double z1 = depth;
        double z0 = std::max(0.0, depth - std::max(0.0, Lm));
        DataConfig cfg = baseCfg_;
        cfg.enh.enable = true; cfg.enh.z_start_m = z0; cfg.enh.z_end_m = z1;
        // Prepare output dir
        std::string outdir = runsRoot + std::string("/") + std::string("L") + std::to_string(int(std::round(Lm))) + std::string("m");
        ensure_dir(outdir);
        std::string resultsPath = outdir + std::string("/results.csv");
        // Run
        SimulationController sim(cfg);
        if (!sim.run(resultsPath)) { _putenv_s("SUMMARY_ONLY","0"); return false; }
        // Summaries
        double comp=0,pump=0,heat=0,LCOH=0; readSummary(outdir+std::string("/summary.csv"), comp,pump,heat,LCOH);
        double elec = comp + pump;
        double scop = elec>1e-9 ? heat/elec : 0.0;
        double unmet = computeUnmetFromResults(resultsPath, dt_h);
        rows.push_back(Row{Lm,z0,z1,comp,pump,heat,elec,scop,LCOH,unmet});
    }

    _putenv_s("SUMMARY_ONLY", "0");
    // Write results
    {
        std::string out = runsRoot + std::string("/opt_results.csv");
        std::ofstream ofs(out.c_str()); if(ofs){
            ofs << "L_enh_m,z_start_m,z_end_m,comp_kWh,pump_kWh,heat_kWh,elec_kWh,SCOP_sys,LCOH_annual,unmet_kWh\n";
            for (auto& r: rows){
                ofs << r.L_enh_m << "," << r.z0 << "," << r.z1 << "," << r.comp << "," << r.pump << "," << r.heat << "," << r.elec << "," << r.scop << "," << r.LCOH << "," << r.unmet << "\n";
            }
        }
    }
    // Pareto (minimize LCOH, elec, unmet)
    {
        std::string outp = runsRoot + std::string("/opt_pareto.csv");
        std::ofstream ofs(outp.c_str()); if(ofs){
            ofs << "L_enh_m,comp_kWh,pump_kWh,heat_kWh,elec_kWh,SCOP_sys,LCOH_annual,unmet_kWh\n";
            for (size_t i=0;i<rows.size();++i){
                auto& a = rows[i]; bool dom=false; for(size_t j=0;j<rows.size();++j){ if(i==j) continue; auto& b = rows[j];
                    bool ge = (b.LCOH <= a.LCOH+1e-9) && (b.elec <= a.elec+1e-9) && (b.unmet <= a.unmet+1e-9);
                    bool gt = (b.LCOH < a.LCOH-1e-9) || (b.elec < a.elec-1e-9) || (b.unmet < a.unmet-1e-9);
                    if (ge && gt){ dom=true; break; }
                }
                if(!dom){ ofs << a.L_enh_m << "," << a.comp << "," << a.pump << "," << a.heat << "," << a.elec << "," << a.scop << "," << a.LCOH << "," << a.unmet << "\n"; }
            }
        }
    }
    return true;
}
