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
    double dt_s,
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
            // Clamp cycle temperatures to CoolProp-safe domain
            // Query fluid critical and triple temperatures to avoid out-of-range calls
            double Tcrit = 0.0, Ttrip = 0.0;
            try {
                Tcrit = g_coolprop.PropsSI("Tcrit", "", 0.0, "", 0.0, fluid.c_str());
                Ttrip = g_coolprop.PropsSI("T_triple", "", 0.0, "", 0.0, fluid.c_str());
            } catch(...) { /* leave zeros, guard below */ }
            if (!(Tcrit > 0.0)) Tcrit = 1e9; // very high sentinel
            if (!(Ttrip > 0.0)) Ttrip = 50.0; // ~50 K sentinel

            // Enforce ordering and margins between evap/cond temperatures
            constexpr double Kmin = 0.5; // minimal approach margin in K
            double Tevap_sat = Tevap_sat_guard_K;
            double Tcond_sat = Tcond_sat_K;
            // Boundaries
            Tevap_sat = std::max(Ttrip + Kmin, std::min(Tevap_sat, Tcrit - 5.0));
            Tcond_sat = std::max(Tevap_sat + 2.0*Kmin, std::min(Tcond_sat, Tcrit - Kmin));

            auto P_evap_Pa = g_coolprop.PropsSI("P","T",Tevap_sat,"Q",1.0,fluid.c_str());
            auto P_cond_Pa = g_coolprop.PropsSI("P","T",Tcond_sat,         "Q",0.0,fluid.c_str());
            last_debug_.P_evap_kPa = P_evap_Pa/1000.0;
            last_debug_.P_cond_kPa = P_cond_Pa/1000.0;
            auto T_superheat_K = std::max(Tevap_sat + Kmin, Tevap_sat + cfg_.hp.superheat_K);
            auto h1 = g_coolprop.PropsSI("H","T",T_superheat_K,"P",P_evap_Pa,fluid.c_str());
            last_debug_.h1 = h1;
            auto s1 = g_coolprop.PropsSI("S","T",T_superheat_K,"P",P_evap_Pa,fluid.c_str());
            auto h2s = g_coolprop.PropsSI("H","P",P_cond_Pa,   "S",s1,       fluid.c_str());
            auto h2  = h1 + (h2s - h1) / (std::max)(0.1, cfg_.hp.eta_isentropic);
            last_debug_.h2s = h2s;
            last_debug_.h2  = h2;
            auto T_subcool_K = (std::max)(Ttrip + Kmin, Tcond_sat - std::max(1.0, cfg_.hp.subcool_K));
            if (T_subcool_K > Tcond_sat - Kmin) T_subcool_K = Tcond_sat - Kmin;
            auto h3 = g_coolprop.PropsSI("H","T",T_subcool_K,"P",P_cond_Pa,fluid.c_str());
            last_debug_.h3 = h3;
            double num = (h2 - h3);
            double den = (std::max)(1e-6, h2 - h1);
            COP = num/den;
            // Sanity bound COP with Carnot-based efficiency ceiling and a hard cap
            double T_evap_K = T_superheat_K;
            double T_cond_K = T_subcool_K + std::max(1.0, cfg_.hp.subcool_K);
            if (T_cond_K - T_evap_K < 1.0) T_evap_K = T_cond_K - 1.0;
            double COP_carnot = (T_cond_K > T_evap_K + 1e-6) ? (T_cond_K / (T_cond_K - T_evap_K)) : 1.0;
            double eff = cfg_.hp.eff_carnot * (0.9 * cfg_.hp.eta_isentropic + 0.1);
            // remove COP upper cap; keep lower bound for stability
            if (!std::isfinite(COP) || COP <= 0.0) COP = 1.0;
            COP = (std::max)(1.0, COP);
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
        // remove COP upper cap in fallback
        COP = (std::max)(1.5, eff * COP_carnot);
    }

    double Q_target = (std::max)(0.0, cfg_.hp.Q_demand_kW);
    const double mod_span = (std::max)(1.0, cfg_.hp.mod_range_C);
    double scale = (T_source_out_C - T_freeze_C) / mod_span;
    scale = clampd(scale, 0.0, 1.0);
    Q_target *= scale;

    // If no demand requested (e.g., compressor off), return no output and zero power
    if (Q_target <= 1e-9) {
        Q_out_kW = 0.0;
        P_el_kW  = 0.0;
        COP      = 0.0;
        return T_source_out_C;
    }

    const double m_dot = (std::max)(1e-9, cfg_.fluid.massFlow_kgps);
    const double cp    = cfg_.fluid.cp;

    // Compute capacity limit from source dT cap and nameplate
    // If HP_MAX_SRC_DT_PER_H <= 0, disable this cap and rely on freeze guard only.
    double dt_hr = (dt_s > 0.0 ? (dt_s / 3600.0) : 0.0);
    double dT_cap = 0.0; // K per step
    if (cfg_.hp.max_source_dT_per_hour > 0.0) {
        dT_cap = cfg_.hp.max_source_dT_per_hour * (dt_hr > 0.0 ? dt_hr : 1.0);
    } else {
        dT_cap = 1e9; // effectively no cap; limited by freeze guard below
    }
    // also enforce freeze guard safety on this step
    double max_drop_to_guard = std::max(0.0, T_source_out_C - cfg_.hp.min_source_return_C);
    dT_cap = std::min(dT_cap, max_drop_to_guard);
    if (COP <= 1.0 || dT_cap <= 0.0) {
        Q_out_kW = 0.0; P_el_kW = 0.0; return T_source_out_C; // no capacity
    }
    const double Q_out_cap_src_kW = (COP / (COP - 1.0)) * (m_dot * cp * dT_cap) / 1000.0;

    // Demand after modulation (freeze-guard scale already applied in Q_target)
    double Q_out_demand_kW = Q_target;
    // Final served output limited by source capacity + nameplate
    Q_out_kW = std::max(0.0, std::min({ Q_out_demand_kW, Q_out_cap_src_kW, cfg_.hp.max_Q_out_kW }));
    P_el_kW  = Q_out_kW / COP;

    // Back-compute source dT actually used
    double Q_geo_served_kW = Q_out_kW * (1.0 - 1.0 / COP);
    double dT_need = (Q_geo_served_kW * 1000.0) / (m_dot * cp);
    dT_need = clampd(dT_need, 0.0, dT_cap);

    double T_return_C = T_source_out_C - dT_need;
    if (T_return_C < cfg_.hp.min_source_return_C) T_return_C = cfg_.hp.min_source_return_C;

    return T_return_C;
}





