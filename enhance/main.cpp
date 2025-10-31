﻿#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#include "DataConfig.h"
#include "SimulationController.h"
#include "Optimizer.h"
#include <cstdlib>
#include <algorithm>
#include <string>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
    // Set default OpenMP threads to 50% of logical processors unless user overrides via OMP_NUM_THREADS
#ifdef _OPENMP
    // 榛樿鍥哄畾 16 绾跨▼锛涘璁剧疆浜嗕互涓嬬幆澧冨彉閲忓垯瑕嗙洊锛?
    //   OMP_THREADS锛堢嚎绋嬫暟锛夋垨 OMP_NUM_THREADS锛堟爣鍑?OpenMP 鍙橀噺锛?
    //   OMP_THREADS_PCT锛堟寜鎬婚€昏緫鏍哥櫨鍒嗘瘮 1..100 璁＄畻锛?
    int threads = 12;
    if (const char* v = std::getenv("OMP_THREADS")) {
        try { int t = std::stoi(v); if (t > 0) threads = t; } catch(...) {}
    } else if (const char* v = std::getenv("OMP_NUM_THREADS")) {
        try { int t = std::stoi(v); if (t > 0) threads = t; } catch(...) {}
    } else if (const char* v = std::getenv("OMP_THREADS_PCT")) {
        try {
            int pct = std::stoi(v);
            if (pct > 0 && pct <= 100) {
                int procs = 1; try { procs = omp_get_num_procs(); } catch(...) { procs = 1; }
                long long t = (static_cast<long long>(procs) * pct + 99) / 100;
                threads = (int)std::max(1LL, t);
            }
        } catch(...) {}
    }
    if (threads < 1) threads = 1;
    omp_set_dynamic(0);            // 鍥哄畾绾跨▼鏁帮紝閬垮厤鑷姩缂╂斁瀵艰嚧鍗犵敤蹇介珮蹇戒綆
    omp_set_num_threads(threads);
#endif
    // Non-interactive defaults; allow overrides via environment variables
    if (const char* v = std::getenv("HP_USE_COOLPROP")) { CFG.hp.use_coolprop = (*v!='0' && *v!='n' && *v!='N'); }
    if (const char* v = std::getenv("HP_FLUID"))        { CFG.hp.fluid = v; }
    if (const char* v = std::getenv("HP_ETA_ISEN"))     try { CFG.hp.eta_isentropic = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_SUPERHEAT"))    try { CFG.hp.superheat_K   = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_SUBCOOL"))      try { CFG.hp.subcool_K    = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_EVAP_APPROACH_K")) try { CFG.hp.evap_approach_K = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_COND_APPROACH_K")) try { CFG.hp.cond_approach_K = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_EFF_CARNOT"))      try { CFG.hp.eff_carnot      = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_MOD_RANGE_C"))     try { CFG.hp.mod_range_C     = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_MAX_Q_OUT_KW")) try { CFG.hp.max_Q_out_kW  = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_MAX_SRC_DT_PER_H")) try { CFG.hp.max_source_dT_per_hour = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_MIN_SRC_RETURN_C")) try { CFG.hp.min_source_return_C    = std::stod(v); } catch(...) {}
    // Default: enable simple weather-driven load using weather.csv in working directory
    CFG.load.enable_weather = true;
    CFG.load.weather_csv = "weather.csv";
    if (const char* v = std::getenv("LOAD_WEATHER_CSV")){ CFG.load.enable_weather=true; CFG.load.weather_csv=v; }
    if (const char* v = std::getenv("LOAD_COLUMN"))     { CFG.load.column_name=v; }
    if (const char* v = std::getenv("LOAD_UA_KW_PER_K"))try { CFG.load.UA_kW_per_K=std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("LOAD_BASE_KW"))    try { CFG.load.base_kW    = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("LOAD_INDOOR_T"))   try { CFG.load.indoor_T_C = std::stod(v);} catch(...) {}

    // Tank configuration overrides
    if (const char* v = std::getenv("TANK_VOL_M3"))     try { CFG.tank.volume_m3   = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("TANK_SET_C"))      try { CFG.tank.setpoint_C  = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("TANK_DB_K"))       try { CFG.tank.deadband_K  = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("TANK_MIN_C"))      try { CFG.tank.min_T_C     = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("TANK_UA_KW_PER_K"))try { CFG.tank.UA_kW_per_K = std::stod(v);} catch(...) {}

    // DHW overrides
    if (const char* v = std::getenv("DHW_ENABLE"))      { CFG.dhw.enable = (*v!='0' && *v!='n' && *v!='N'); }
    if (const char* v = std::getenv("DHW_BASE_KW"))     try { CFG.dhw.base_kW      = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("DHW_MORN_START"))  try { CFG.dhw.morning_start_h = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("DHW_MORN_HOURS"))  try { CFG.dhw.morning_hours   = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("DHW_MORN_KW"))     try { CFG.dhw.morning_kW      = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("DHW_EVE_START"))   try { CFG.dhw.evening_start_h = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("DHW_EVE_HOURS"))   try { CFG.dhw.evening_hours   = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("DHW_EVE_KW"))      try { CFG.dhw.evening_kW      = std::stod(v);} catch(...) {}

    // Pump hydraulics overrides
    if (const char* v = std::getenv("PUMP_EFF"))           try { CFG.pump.eff = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("PUMP_REL_ROUGH"))     try { CFG.pump.rel_roughness = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("PUMP_K_MINOR_ANN"))   try { CFG.pump.K_minor_annulus = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("PUMP_K_MINOR_INNER")) try { CFG.pump.K_minor_inner   = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("MASS_FLOW_KGPS"))     try { CFG.fluid.massFlow_kgps  = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("H_IO_WM2K"))          try { CFG.h_io_Wm2K            = std::stod(v);} catch(...) {}

    // Pipe material: set PE-like conductivity by default (reduce short-circuiting)
    // Thickness stays unchanged; override via PIPE_K if needed.
    CFG.pipe.k = 0.4; // W/m-K (typical PE/PEX ~0.35鈥?.45)
    if (const char* v = std::getenv("PIPE_K"))             try { CFG.pipe.k = std::stod(v);} catch(...) {}

    // Geometry and ground overrides (optional)
    // Geometry and ground overrides (optional)
    // Inner/outer pipe conductivity overrides
    if (const char* v = std::getenv("PIPE_K_INNER"))    try { CFG.pipe_k_inner = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("PIPE_K_OUTER"))    try { CFG.pipe_k_outer = std::stod(v);} catch(...) {}
    // Inner-pipe top insulation
    if (const char* v = std::getenv("INSUL_TOP_LEN_M")) try { CFG.insul.top_len_m = std::stod(v); CFG.insul.enable = (CFG.insul.top_len_m>0.0);} catch(...) {}
    if (const char* v = std::getenv("INSUL_K_INNER"))   try { CFG.insul.k_inner   = std::stod(v);} catch(...) {}
    // Bottom-only enhancement length (from bottom up)
    if (const char* v = std::getenv("EHEP_BOTTOM_LEN_M")) { try { double len=std::stod(v); if (len>0) { double L=CFG.well.depth_m; CFG.enh.enable=true; CFG.enh.z_start_m=std::max(0.0, L-len); CFG.enh.z_end_m=L; } } catch(...) {} }
    if (const char* v = std::getenv("WELL_DEPTH_M"))       try { CFG.well.depth_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("WELL_SEGS"))          try { CFG.well.segments = std::max(1, std::stoi(v)); } catch(...) {}
    if (const char* v = std::getenv("D_OUTER_M"))          try { CFG.well.D_outer_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("D_INNER_M"))          try { CFG.well.D_inner_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("BORE_D_M"))           try { CFG.well.borehole_D_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("PIPE_THICK_M"))       try { CFG.well.pipe_wall_thickness_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("T_SURF_C"))           try { CFG.T_surface_C = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("GEO_GRAD_CPM"))       try { CFG.geograd_C_per_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("SOIL_CP"))            try { CFG.soil.cp = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("SOIL_K"))             try { CFG.soil.k = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("SOIL_RHO"))           try { CFG.soil.rho = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("GROUT_K"))            try { CFG.grout.k = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("GROUT_RHO"))          try { CFG.grout.rho = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("GROUT_CP"))           try { CFG.grout.cp = std::stod(v);} catch(...) {}

    // Enhanced pipe overrides
    if (const char* v = std::getenv("EHEP_ENABLE"))     { CFG.enh.enable = (*v!='0' && *v!='n' && *v!='N'); }
    if (const char* v = std::getenv("EHEP_NU_MULT"))    try { CFG.enh.nu_mult = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("EHEP_F_MULT"))     try { CFG.enh.f_mult  = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("EHEP_Z_START_M")) try { CFG.enh.z_start_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("EHEP_Z_END_M"))   try { CFG.enh.z_end_m   = std::stod(v);} catch(...) {}

    // Economics overrides
    if (const char* v = std::getenv("ECON_ELEC_PRICE"))   try { CFG.econ.elec_price_per_kWh = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_CAPEX"))        try { CFG.econ.capex = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_LIFETIME_Y"))   try { CFG.econ.lifetime_years = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_DR"))           try { CFG.econ.discount_rate = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_OM_FACTOR"))    try { CFG.econ.om_factor = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_DRILL_COST_PER_M")) try { CFG.econ.drill_cost_per_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_CASING_COST_PER_M")) try { CFG.econ.casing_cost_per_m = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_EHEP_COST_PER_M"))   try { CFG.econ.ehep_cost_per_m   = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_HP_UNIT_COST_PER_KW"))   try { CFG.econ.hp_unit_cost_per_kW   = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_PUMP_UNIT_COST_PER_KW")) try { CFG.econ.pump_unit_cost_per_kW = std::stod(v);} catch(...) {}

    // Calendar & heating season overrides
    if (const char* v = std::getenv("SIM_START_YEAR"))   try { CFG.sim_start_year  = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("SIM_START_MONTH"))  try { CFG.sim_start_month = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("SIM_START_DAY"))    try { CFG.sim_start_day   = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("HEAT_SEASON_ENABLE")) { CFG.season.enable = (*v!='0' && *v!='n' && *v!='N'); }
    if (const char* v = std::getenv("HEAT_START_MM"))    try { CFG.season.start_month = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("HEAT_START_DD"))    try { CFG.season.start_day   = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("HEAT_END_MM"))      try { CFG.season.end_month   = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("HEAT_END_DD"))      try { CFG.season.end_day     = std::stoi(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_CASING_COST_PER_M"))try { CFG.econ.casing_cost_per_m= std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_EHEP_COST_PER_M"))  try { CFG.econ.ehep_cost_per_m  = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_HP_COST"))         try { CFG.econ.hp_cost_fixed     = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_TANK_COST"))       try { CFG.econ.tank_cost_fixed   = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_PUMP_COST"))       try { CFG.econ.pump_cost_fixed   = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("ECON_MISC_FIXED"))      try { CFG.econ.misc_fixed        = std::stod(v);} catch(...) {}

    // Simulation horizon: default to 1 year at 1-hour timestep unless overridden
    {
        int years = 1;
        if (const char* v = std::getenv("SIM_YEARS")) { try { years = std::max(1, std::stoi(v)); } catch(...) {} }
        int dt_s = 600; // default 10 minutes per step
        if (const char* v = std::getenv("TIME_STEP_S")) { try { dt_s = std::max(60, std::stoi(v)); } catch(...) {} }
        CFG.time.timeStep_s = static_cast<double>(dt_s);
        // steps per hour = 3600/dt_s
        int steps_per_year = static_cast<int>( (3600.0 / CFG.time.timeStep_s) * 8760.0 + 0.5 );
        CFG.time.totalSteps = years * steps_per_year;
    }

    // Default bottom-only enhancement zone if enabled but no interval provided
    if (CFG.enh.enable) {
        if (!(CFG.enh.z_end_m > CFG.enh.z_start_m)) {
            double L = CFG.well.depth_m;
            CFG.enh.z_start_m = 0.7 * L; // bottom 30% by default
            CFG.enh.z_end_m   = L;
        }
    }

    // Multi-segment enhancement: optional CSV, else fallback to single segment
        if (const char* path = std::getenv("EHEP_SEGMENTS_CSV")){
            FILE* f = fopen(path, "r");
            if (f){ char buf[512]; bool first=true; while(fgets(buf,sizeof(buf),f)){
                std::string line(buf); if(first){ first=false; if(line.find("z_start")!=std::string::npos) continue; }
                std::stringstream ss(line); std::string c; std::vector<double> v; while(std::getline(ss,c,',')){ try{ v.push_back(std::stod(c)); } catch(...){} }
                if (v.size()>=4){ EHEPSegment s; s.z_start_m=v[0]; s.z_end_m=v[1]; s.nu_mult=v[2]; s.f_mult=v[3]; CFG.enh_segments.push_back(s);} }
                fclose(f);
            }
        }
    {
        auto baseNu = CFG.NuFunc;
        CFG.NuFunc = [baseNu](double Re, double Pr, double z)->double{
            double nu = 50.0; try { if (baseNu) nu = baseNu(Re, Pr, z); } catch(...) {}
            double mult = 1.0;
            if (!CFG.enh_segments.empty()){
                for (auto& s: CFG.enh_segments){ if (z>=s.z_start_m && z<=s.z_end_m) mult *= std::max(1.0, s.nu_mult); }
            } else if (CFG.enh.enable){
                bool enhAll = (CFG.enh.z_end_m <= CFG.enh.z_start_m) || (CFG.enh.z_end_m <= 0.0);
                if (enhAll || (z >= CFG.enh.z_start_m && z <= CFG.enh.z_end_m)) mult *= std::max(1.0, CFG.enh.nu_mult);
            }
            return nu * mult;
        };
    }

    // Integrated optimizer switch
    const char* runopt = std::getenv("RUN_OPTIMIZE");
    if (runopt && (*runopt!='0' && *runopt!='n' && *runopt!='N')){
        // Runs root dir
        std::string runsRoot = "opt_runs";
        if (const char* v = std::getenv("OPT_DIR")) runsRoot = v;
        Optimizer opt(CFG);
        bool ok = opt.run(runsRoot);
        if (!ok) return 2;
    } else {
        SimulationController sim(CFG);
        sim.run("results.csv");
    }
    return 0;
}


