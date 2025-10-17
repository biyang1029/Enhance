#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#include "DataConfig.h"
#include "SimulationController.h"
#include <cstdlib>\n#include <string>

int main() {
    // Non-interactive defaults; allow overrides via environment variables
    if (const char* v = std::getenv("HP_USE_COOLPROP")) { CFG.hp.use_coolprop = (*v!='0' && *v!='n' && *v!='N'); }
    if (const char* v = std::getenv("HP_FLUID"))        { CFG.hp.fluid = v; }
    if (const char* v = std::getenv("HP_ETA_ISEN"))     try { CFG.hp.eta_isentropic = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_SUPERHEAT"))    try { CFG.hp.superheat_K   = std::stod(v); } catch(...) {}
    if (const char* v = std::getenv("HP_SUBCOOL"))      try { CFG.hp.subcool_K    = std::stod(v); } catch(...) {}
    // Default: enable simple weather-driven load using weather.csv in working directory
    CFG.load.enable_weather = true;
    CFG.load.weather_csv = "weather.csv";
    if (const char* v = std::getenv("LOAD_WEATHER_CSV")){ CFG.load.enable_weather=true; CFG.load.weather_csv=v; }
    if (const char* v = std::getenv("LOAD_COLUMN"))     { CFG.load.column_name=v; }
    if (const char* v = std::getenv("LOAD_UA_KW_PER_K"))try { CFG.load.UA_kW_per_K=std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("LOAD_BASE_KW"))    try { CFG.load.base_kW    = std::stod(v);} catch(...) {}
    if (const char* v = std::getenv("LOAD_INDOOR_T"))   try { CFG.load.indoor_T_C = std::stod(v);} catch(...) {}

    SimulationController sim(CFG);
    sim.run("results.csv");
    return 0;
}

