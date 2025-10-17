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

    SimulationController sim(CFG);
    sim.run("results.csv");
    return 0;
}

