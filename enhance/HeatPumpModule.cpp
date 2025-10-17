#include "HeatPumpModule.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#include <windows.h>
#endif

static inline double clampd(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

#ifdef _WIN32
struct CoolPropDyn {
    HMODULE h = nullptr;
    typedef double (__cdecl *PropsSI_t)(const char*, const char*, double, const char*, double, const char*);
    PropsSI_t PropsSI = nullptr;
    bool load_once = false;
    bool try_load(const std::string& hintDir = std::string()){
        if (load_once) return h && PropsSI;
        load_once = true;
        auto try_path = [&](const std::string& path){
            h = LoadLibraryA(path.c_str());
            if (h){ PropsSI = (PropsSI_t)GetProcAddress(h, "PropsSI"); }
            return h && PropsSI;
        };
        // 1) explicit hint
        if (!hintDir.empty()){
            std::string p = hintDir;
            char back = p.empty()? '\0' : p.back();
            if (back!='\\' && back!='/') p += "\\";
            if (try_path(p + "CoolProp.dll")) return true;
        }
        // 2) COOLPROP_DIR env (bin and root)
        char* env = nullptr; size_t len=0;
        if (_dupenv_s(&env,&len,"COOLPROP_DIR")==0 && env){
            std::string base(env); free(env);
            std::string p1 = base + "\\bin\\CoolProp.dll";
            if (try_path(p1)) return true;
            std::string p2 = base + "\\CoolProp.dll";
            if (try_path(p2)) return true;
        }
        // 3) current dir
        if (try_path("CoolProp.dll")) return true;
        // 4) common python wheel location (Anaconda)
        if (_dupenv_s(&env,&len,"USERPROFILE")==0 && env){
            std::string up(env); free(env);
            // Best-effort guess; user can still drop DLL next to exe
            std::string guess = "D:/yangb/soft/anaconda/Lib/site-packages/CoolProp/CoolProp.dll";
            try_path(guess);
        }
        return h && PropsSI;
    }
};
static CoolPropDyn g_coolprop;
#endif

double HeatPumpModule::step(double T_source_out_C,
    double& COP,
    double& Q_out_kW,
    double& P_el_kW) {
    const double T_load_set_C = cfg_.hp.target_T_load_out_C;

    const double Tevap_sat_K = (T_source_out_C - cfg_.hp.evap_approach_K) + 273.15;
    const double Tcond_sat_K = (T_load_set_C   + cfg_.hp.cond_approach_K) + 273.15;

    const double T_freeze_C = cfg_.hp.freeze_guard_C;
    const double Tevap_sat_guard_K = (std::max)(Tevap_sat_K, (T_freeze_C + 0.1) + 273.15);

    bool used_cp = false;
    // reset debug
    last_debug_ = HPDebug{};
    last_debug_.fluid = cfg_.hp.fluid;
    last_debug_.T_evap_sat_K = Tevap_sat_K;
    last_debug_.T_cond_sat_K = Tcond_sat_K;
    static bool warned_failure = false;

#ifdef _WIN32
    if (cfg_.hp.use_coolprop && g_coolprop.try_load()){
        try {
            const std::string& fluid = cfg_.hp.fluid;
            auto P_evap_Pa = g_coolprop.PropsSI("P","T",Tevap_sat_guard_K,"Q",1.0,fluid.c_str());
            auto P_cond_Pa = g_coolprop.PropsSI("P","T",Tcond_sat_K,         "Q",0.0,fluid.c_str());
            last_debug_.P_evap_kPa = P_evap_Pa/1000.0;
            last_debug_.P_cond_kPa = P_cond_Pa/1000.0;
            auto T_superheat_K = Tevap_sat_guard_K + cfg_.hp.superheat_K;
            auto h1 = g_coolprop.PropsSI("H","T",T_superheat_K,"P",P_evap_Pa,fluid.c_str());
            last_debug_.h1 = h1;
            auto s1 = g_coolprop.PropsSI("S","T",T_superheat_K,"P",P_evap_Pa,fluid.c_str());
            auto h2s = g_coolprop.PropsSI("H","P",P_cond_Pa,   "S",s1,       fluid.c_str());
            auto h2  = h1 + (h2s - h1) / (std::max)(0.1, cfg_.hp.eta_isentropic);
            last_debug_.h2s = h2s;
            last_debug_.h2  = h2;
            auto T_subcool_K = (std::max)(1.0, Tcond_sat_K - cfg_.hp.subcool_K);
            auto h3 = g_coolprop.PropsSI("H","T",T_subcool_K,"P",P_cond_Pa,fluid.c_str());
            last_debug_.h3 = h3;
            double num = (h2 - h3);
            double den = (std::max)(1e-6, h2 - h1);
            COP = clampd(num/den, 1.0, 10.0);
            used_cp = std::isfinite(COP);
            last_debug_.used_coolprop = used_cp;
            if (!used_cp && !warned_failure){
                std::cerr << "[Warning] CoolProp PropsSI returned non-finite values. Fallback model used.\n";
                warned_failure = true;
            }
        } catch(...) {
            if (!warned_failure){ std::cerr << "[Warning] CoolProp call threw. Fallback model used.\n"; warned_failure = true; }
        }
    }
#endif

    if (!used_cp){
        double T_evap_K = Tevap_sat_guard_K + cfg_.hp.superheat_K;
        double T_cond_K = Tcond_sat_K - cfg_.hp.subcool_K;
        const double dT_min = 5.0;
        if (T_cond_K - T_evap_K < dT_min) T_evap_K = T_cond_K - dT_min;
        double COP_carnot = (T_cond_K > T_evap_K + 1e-6) ? (T_cond_K / (T_cond_K - T_evap_K)) : 1.0;
        double eff = cfg_.hp.eff_carnot * (0.9 * cfg_.hp.eta_isentropic + 0.1);
        COP = clampd(eff * COP_carnot, 1.5, 8.0);
    }

    double Q_target = (std::max)(0.0, cfg_.hp.Q_demand_kW);
    const double mod_span = (std::max)(1.0, cfg_.hp.mod_range_C);
    double scale = (T_source_out_C - T_freeze_C) / mod_span;
    scale = clampd(scale, 0.0, 1.0);
    Q_target *= scale;

    const double m_dot = (std::max)(1e-9, cfg_.fluid.massFlow_kgps);
    const double cp    = cfg_.fluid.cp;

    const double Q_geo_need_kW = (COP > 1.0) ? Q_target * (1.0 - 1.0 / COP) : 0.0;
    double dT_need = (Q_geo_need_kW * 1000.0) / (m_dot * cp);

    const double dT_max = (std::max)(0.5, cfg_.hp.max_source_dT_per_hour);
    dT_need = clampd(dT_need, 0.0, dT_max);

    double T_return_C = T_source_out_C - dT_need;
    if (T_return_C < T_freeze_C) {
        T_return_C = T_freeze_C;
        dT_need = T_source_out_C - T_return_C;
    }

    const double Q_geo_act_kW = (m_dot * cp * dT_need) / 1000.0;
    Q_out_kW = (COP > 1.0) ? (COP / (COP - 1.0)) * Q_geo_act_kW : 0.0;
    P_el_kW  = (COP > 0.0) ? Q_out_kW / COP : 0.0;

    return T_return_C;
}

