#include "HeatPumpModule.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include "RefpropWrapper.h"

static inline double clampd(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

double HeatPumpModule::step(double T_source_out_C,
    double& COP,
    double& Q_out_kW,
    double& P_el_kW) {
    const double T_load_set_C = cfg_.hp.target_T_load_out_C;

    // State-point temperatures (sat + approach + superheat/subcool)
    const double Tevap_sat_K = (T_source_out_C - cfg_.hp.evap_approach_K) + 273.15;
    const double Tcond_sat_K = (T_load_set_C   + cfg_.hp.cond_approach_K) + 273.15;

    // Guard freezing
    const double T_freeze_C = cfg_.hp.freeze_guard_C;
    const double Tevap_sat_guard_K = std::max(Tevap_sat_K, (T_freeze_C + 0.1) + 273.15);

    // Try REFPROP first (optional)
    bool used_refprop = false;
    if (cfg_.hp.use_refprop) {
        static RefpropWrapper rp;
        static bool load_attempted = false;
        static bool load_success = false;
        static std::string last_dir;
        static bool warned_failure = false;

        // Reload if directory changed
        if (!load_attempted || cfg_.hp.refprop_dll_dir != last_dir) {
            last_dir = cfg_.hp.refprop_dll_dir;
            load_success = rp.load(last_dir);
            load_attempted = true;
            if (!load_success && !warned_failure) {
                std::cerr << "[Warning] REFPROP load failed from \"" << last_dir
                          << "\". Falling back to empirical COP model.\n";
                warned_failure = true;
            }
        }

        if (load_success) {
            rp.setFluid(cfg_.hp.fluid);
            double P_evap_kPa = 0.0;
            double P_cond_kPa = 0.0;
            bool ok = rp.psat_T(Tevap_sat_guard_K, P_evap_kPa, 1) &&
                      rp.psat_T(Tcond_sat_K, P_cond_kPa, 1);
            if (ok) {
                double h1 = 0.0;
                double s1 = 0.0;
                double h2s = 0.0;
                double T2s = 0.0;
                double h3 = 0.0;
                double s3 = 0.0;

                ok = rp.hs_from_TP(Tevap_sat_guard_K + cfg_.hp.superheat_K,
                                   P_evap_kPa, h1, s1);
                ok = ok && rp.hT_from_PS(P_cond_kPa, s1, h2s, T2s);
                ok = ok && rp.hs_from_TP(std::max(1.0,
                                   Tcond_sat_K - cfg_.hp.subcool_K),
                                   P_cond_kPa, h3, s3);

                if (ok) {
                    double h2 = h1 + (h2s - h1) / std::max(0.1, cfg_.hp.eta_isentropic);
                    double den = std::max(1e-6, h2 - h1);
                    double num = (h2 - h3);
                    COP = clampd(num / den, 1.0, 10.0);
                    used_refprop = true;
                } else if (!warned_failure) {
                    std::cerr << "[Warning] REFPROP evaluation failed for fluid " << cfg_.hp.fluid
                              << ". Falling back to empirical COP model.\n";
                    warned_failure = true;
                }
            } else if (!warned_failure) {
                std::cerr << "[Warning] REFPROP saturation call failed for fluid " << cfg_.hp.fluid
                          << ". Falling back to empirical COP model.\n";
                warned_failure = true;
            }
        }
    }

    if (!used_refprop){
        // fallback: approach + Carnot scaling + eta_isentropic factor
        double T_evap_K = Tevap_sat_guard_K + cfg_.hp.superheat_K;
        double T_cond_K = Tcond_sat_K - cfg_.hp.subcool_K;
        const double dT_min = 5.0;
        if (T_cond_K - T_evap_K < dT_min) T_evap_K = T_cond_K - dT_min;
        double COP_carnot = (T_cond_K > T_evap_K + 1e-6) ? (T_cond_K / (T_cond_K - T_evap_K)) : 1.0;
        double eff = cfg_.hp.eff_carnot * (0.9 * cfg_.hp.eta_isentropic + 0.1);
        COP = clampd(eff * COP_carnot, 1.5, 8.0);
    }

    // Target load with freeze guard modulation
    double Q_target = std::max(0.0, cfg_.hp.Q_demand_kW);
    const double mod_span = std::max(1.0, cfg_.hp.mod_range_C);
    double scale = (T_source_out_C - T_freeze_C) / mod_span;
    scale = clampd(scale, 0.0, 1.0);
    Q_target *= scale;

    // Source-side deltaT constraint
    const double m_dot = std::max(1e-9, cfg_.fluid.massFlow_kgps);
    const double cp    = cfg_.fluid.cp;

    const double Q_geo_need_kW = (COP > 1.0) ? Q_target * (1.0 - 1.0 / COP) : 0.0;
    double dT_need = (Q_geo_need_kW * 1000.0) / (m_dot * cp);

    const double dT_max = std::max(0.5, cfg_.hp.max_source_dT_per_hour);
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
