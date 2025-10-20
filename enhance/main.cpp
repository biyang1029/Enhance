#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#include "DataConfig.h"
#include "SimulationController.h"
#include <cstdlib>
#include <string>

int main() {
    // Non-interactive defaults; allow overrides via environment variables
    if (const char* v = std::getenv("HP_USE_COOLPROP")) { CFG.hp.use_coolprop = (*v!='0' && *v!='n' && *v!='N'); }
    if (const char* v = std::getenv("HP_FLUID"))        { CFG.hp.fluid = v; }
    if (const char* v = std::getenv("HP_ETA_ISEN"))     try { CFG.hp.eta_isentropic = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_SUPERHEAT"))    try { CFG.hp.superheat_K   = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_SUBCOOL"))      try { CFG.hp.subcool_K    = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_MAX_Q_OUT_KW")) try { CFG.hp.max_Q_out_kW  = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_MAX_SRC_DT_PER_H")) try { CFG.hp.max_source_dT_per_hour = std::stod(v); } catch(...) {}
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

    // Simulation horizon: default to 1 year at 1-hour timestep unless overridden
    {
        int years = 1;
        if (const char* v = std::getenv("SIM_YEARS")) { try { years = std::max(1, std::stoi(v)); } catch(...) {} }
        int dt_s = 3600;
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

    // Wrap Nu with enhancement multiplier if enabled
    if (CFG.enh.enable) {
        auto baseNu = CFG.NuFunc;
        CFG.NuFunc = [baseNu](double Re, double Pr, double z)->double{
            double nu = 50.0;
            try {
                if (baseNu) nu = baseNu(Re, Pr, z);
            } catch(...) {}
            bool enhAll = (CFG.enh.z_end_m <= CFG.enh.z_start_m) || (CFG.enh.z_end_m <= 0.0);
            if (CFG.enh.enable && (enhAll || (z >= CFG.enh.z_start_m && z <= CFG.enh.z_end_m))) {
                nu *= std::max(1.0, CFG.enh.nu_mult);
            }
            return nu;
        };
    }

    SimulationController sim(CFG);
    sim.run("results.csv");
    return 0;
}

