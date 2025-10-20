// ========================= SimulationController.cpp =========================
#define _CRT_SECURE_NO_WARNINGS 1
#include "SimulationController.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <fstream>

// 进度条打印
static void printProgress(int done, int total, int barWidth = 40) {
    if (done < 0) done = 0;
    if (done > total) done = total;
    double frac = static_cast<double>(done) / static_cast<double>(total);
    int filled = static_cast<int>(frac * barWidth + 0.5);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) std::cout << (i < filled ? '=' : ' ');
    std::cout << "] " << int(frac * 100.0) << "% (" << done << "/" << total << ")" << std::flush;
}

bool SimulationController::run(const std::string& csv_path) {
    // 初始化模块
    ground_.initialize(cfg_);
    hp_.initialize(cfg_);
    tank_.initialize(cfg_);

    if (!logger_.open(csv_path)) {
        std::cerr << "[Error] Failed to open log file: " << csv_path << "\n";
        return false;
    }
    // 同时打开 debug 日志
    logger_.openDebug("debug.csv");

    // 回水初值
    double T_return = cfg_.fluid.inlet_T_C;
    tank_.reset(cfg_.tank.setpoint_C);
    const int    maxIter = 300;
    const double tolT    = 0.005;
    const double relax   = 0.5;

    // Variable flow configuration (min/max kg/s), default around base mass flow
    static bool flow_cfg_inited = false;
    static double flow_min_kgps = 0.0, flow_max_kgps = 0.0;
    if (!flow_cfg_inited) {
        double base = cfg_.fluid.massFlow_kgps;
        flow_min_kgps = std::max(0.1, 0.4 * base);
        flow_max_kgps = std::max(flow_min_kgps, 1.2 * base);
        if (const char* v = std::getenv("FLOW_MIN_KGPS")) { try { flow_min_kgps = std::stod(v); } catch(...) {} }
        if (const char* v = std::getenv("FLOW_MAX_KGPS")) { try { flow_max_kgps = std::stod(v); } catch(...) {} }
        flow_cfg_inited = true;
    }

    // Physical temperature guardrails (avoid runaway states)
    const double T_min_phys = std::min(0.0, cfg_.hp.min_source_return_C - 2.0);
    const double T_max_phys = cfg_.T_surface_C + cfg_.geograd_C_per_m * cfg_.well.depth_m + 60.0;

    for (int h = 0; h < cfg_.time.totalSteps; ++h) {
        double Q_space_req_kW = 0.0;
        // 每小时（如启用）计算天气负荷并注入
        if (cfg_.load.enable_weather) {
            static bool loaded = false;
            if (!loaded) { load_.loadCSV(cfg_.load.weather_csv, cfg_.load.column_name); loaded = true; }
            Q_space_req_kW = load_.demandAt(h, cfg_.load);
        }
        double Q_dhw_req_kW = cfg_.dhw.enable ? load_.dhwAt(h, cfg_.dhw) : 0.0;

        // Tank thermostat on/off with deadband
        const double Ttank = tank_.temperature_C();
        const double on_th  = cfg_.tank.setpoint_C - cfg_.tank.deadband_K;
        const double off_th = cfg_.tank.setpoint_C + cfg_.tank.deadband_K;
        if (!hp_on_ && Ttank <= on_th) hp_on_ = true;
        if (hp_on_  && Ttank >= off_th) hp_on_ = false;

        // Compressor demand based on thermostat
        if (hp_on_) hp_.setDemand(cfg_.hp.max_Q_out_kW); else hp_.setDemand(0.0);

        // Variable flow control: scale between min/max by (space+dhw)/Qmax when ON
        double frac = 0.0;
        if (hp_on_) {
            double qreq = std::max(0.0, Q_space_req_kW + Q_dhw_req_kW);
            double qmax = std::max(1e-6, cfg_.hp.max_Q_out_kW);
            frac = std::clamp(qreq / qmax, 0.0, 1.0);
        }
        double flow_now = flow_min_kgps + (flow_max_kgps - flow_min_kgps) * frac;
        ground_.setMassFlow(flow_now);
        hp_.setMassFlow(flow_now);

        double COP = 0, Q_out_kW = 0, P_el_kW = 0, Q_geo_kW = 0;
        double T_iter = std::clamp(T_return, T_min_phys, T_max_phys);
        double T_source_out = 0.0, T_next = T_iter;

        // 固定点预迭代（不推进土壤）
        for (int it = 0; it < maxIter - 1; ++it) {
            T_source_out = ground_.step(T_iter, Q_geo_kW, /*advanceSoil=*/false);
            if (!std::isfinite(T_source_out)) T_source_out = T_iter;
            T_source_out = std::clamp(T_source_out, T_min_phys, T_max_phys);

            T_next = hp_.step(T_source_out, COP, Q_out_kW, P_el_kW);
            if (!std::isfinite(T_next)) T_next = T_iter;
            // limit per-iteration jump to avoid divergence
            const double max_dT_iter = 10.0; // K per inner iteration
            double dT = std::clamp(T_next - T_iter, -max_dT_iter, max_dT_iter);
            T_next = T_iter + dT;
            T_next = std::clamp(T_next, T_min_phys, T_max_phys);
            if (std::fabs(T_next - T_iter) < tolT) break;
            T_iter = relax * T_next + (1.0 - relax) * T_iter;
        }

        // 最后一次：推进土壤并完成一步
        T_source_out = ground_.step(T_iter, Q_geo_kW, /*advanceSoil=*/true);
        if (!std::isfinite(T_source_out)) T_source_out = T_iter;
        T_source_out = std::clamp(T_source_out, T_min_phys, T_max_phys);
        double T_return_next = hp_.step(T_source_out, COP, Q_out_kW, P_el_kW);
        if (!std::isfinite(T_return_next)) T_return_next = T_iter;
        T_return_next = std::clamp(T_return_next, T_min_phys, T_max_phys);

        // Apply tank energy balance for this step
        double Q_space_served_kW = 0.0, Q_dhw_served_kW = 0.0, Q_unmet_kW = 0.0;
        tank_.applyHour(cfg_.time.timeStep_s, Q_out_kW, Q_space_req_kW, Q_dhw_req_kW,
                        Q_space_served_kW, Q_dhw_served_kW, Q_unmet_kW);

        // Pump hydraulics: Darcy–Weisbach for annulus (down) and inner pipe (up)
        auto safe = [](double v, double lo){ return (v>lo? v: lo); };
        const double rho = cfg_.fluid.rho;
        const double mu  = cfg_.fluid.mu;
        const double L   = std::max(0.0, cfg_.well.depth_m);
        const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
        const double D_o_out = std::max(1e-6, cfg_.well.D_outer_m);
        const double D_o_inn = std::max(1e-6, cfg_.well.D_inner_m);
        const double D_o_inner_wall = std::max(1e-6, D_o_out - 2.0 * t_p); // inner diameter of outer pipe (annulus outer)
        const double D_i_inner = std::max(1e-6, D_o_inn - 2.0 * t_p);      // inner pipe inner diameter (assumed same thickness)
        const double Dh_ann = std::max(1e-6, D_o_inn > D_o_inner_wall ? (D_o_inn - D_o_inner_wall) : (D_o_inner_wall - D_o_inn));
        const double A_ann = std::max(1e-12, 0.25 * PI * (D_o_inner_wall*D_o_inner_wall - D_o_inn*D_o_inn));
        const double A_in  = std::max(1e-12, 0.25 * PI * (D_i_inner*D_i_inner));
        const double m_dot = std::max(1e-9, flow_now);
        const double v_ann = m_dot / (rho * A_ann);
        const double v_in  = m_dot / (rho * A_in);
        const double Re_ann = rho * v_ann * Dh_ann / safe(mu,1e-9);
        const double Re_in  = rho * v_in  * std::max(1e-6, D_i_inner) / safe(mu,1e-9);
        auto fric = [&](double Re, double rel_rough)->double{
            if (Re < 1e-12) return 0.0;
            if (Re < 2300.0) return 64.0 / safe(Re,1e-9);
            // Swamee-Jain explicit
            double term = rel_rough/3.7 + 5.74/std::pow(Re,0.9);
            double f = 0.25/std::pow(std::log10(term),2.0);
            if (!std::isfinite(f) || f<=0.0) f = 0.02; // fallback
            return f;
        };
        const double f_ann = fric(Re_ann, cfg_.pump.rel_roughness);
        const double f_in  = fric(Re_in,  cfg_.pump.rel_roughness);
        const double q_ann = 0.5 * rho * v_ann * v_ann;
        const double q_in  = 0.5 * rho * v_in  * v_in;
        // Apply EHEP friction multiplier on annulus over enhanced fraction of length
        double f_ann_eff = f_ann;
        if (cfg_.enh.enable) {
            double frac_enh = 1.0;
            if (cfg_.enh.z_end_m > cfg_.enh.z_start_m && L > 0.0) {
                frac_enh = std::clamp((cfg_.enh.z_end_m - cfg_.enh.z_start_m) / L, 0.0, 1.0);
            }
            f_ann_eff = f_ann * (1.0 + frac_enh * (std::max(1.0, cfg_.enh.f_mult) - 1.0));
        }
        const double dP_ann = f_ann_eff * (L / std::max(1e-6, Dh_ann)) * q_ann + cfg_.pump.K_minor_annulus * q_ann;
        const double dP_in  = f_in  * (L / std::max(1e-6, D_i_inner)) * q_in  + cfg_.pump.K_minor_inner   * q_in;
        const double dP_total_Pa = dP_ann + dP_in;
        const double eta_p = std::max(0.05, std::min(0.95, cfg_.pump.eff));
        const double P_pump_W = (m_dot / rho) * dP_total_Pa / eta_p;
        const double P_pump_kW = P_pump_W / 1000.0;
        const double dP_kPa = dP_total_Pa / 1000.0;

        // 一次性在控制台打印后端信息
        {
            static bool printed_backend = false;
            auto dbg = hp_.lastDebug();
            if (!printed_backend) {
                std::cout << "\n[Info] Thermo backend: "
                          << (dbg.used_coolprop ? "CoolProp" : "fallback")
                          << " (fluid=" << cfg_.hp.fluid << ")\n";
                printed_backend = true;
            }
        }

        // 写结果 + 模型来源
        {
            auto dbg = hp_.lastDebug();
            std::string model = dbg.used_coolprop ? "CP" : "FB";
            logger_.writeHour(h, T_source_out, T_return_next, COP, Q_out_kW, P_el_kW, Q_geo_kW, model,
                              tank_.temperature_C(), hp_on_?1:0,
                              Q_space_req_kW, Q_dhw_req_kW,
                              Q_space_served_kW, Q_dhw_served_kW, Q_unmet_kW,
                              flow_now, dP_kPa, P_pump_kW);
            logger_.writeDebugHour(h, dbg.fluid, dbg.used_coolprop,
                                   dbg.T_evap_sat_K, dbg.T_cond_sat_K,
                                   dbg.P_evap_kPa, dbg.P_cond_kPa,
                                   dbg.h1, dbg.h2s, dbg.h2, dbg.h3);
        }

        T_return = T_return_next;
        printProgress(h + 1, cfg_.time.totalSteps);
    }

    logger_.close();

    // Economic summary (optional scaling to annual)
    // Re-parse results.csv to aggregate energy if not tracked; here we aggregate from runtime variables
    // For now, recompute by reading results.csv to avoid adding many state vars
    try {
        std::ifstream ifs("results.csv");
        std::string line; if (!ifs) throw 1;
        // header
        std::getline(ifs, line);
        double sumPel_kWh = 0.0, sumPump_kWh = 0.0, sumQspace_kWh = 0.0, sumQdhw_kWh = 0.0;
        const double dt_h = cfg_.time.timeStep_s / 3600.0;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            // crude CSV parse
            // columns: step, T_source_out_C, T_return_C, COP, Q_out_kW, P_el_kW, Q_geo_kW, model, T_tank_C, HP_on, Q_space_req_kW, Q_dhw_req_kW, Q_space_served_kW, Q_dhw_served_kW, Q_unmet_kW, flow_kgps, dP_kPa, P_pump_kW
            // indexes: 0        1                2            3    4         5          6        7      8          9      10               11               12                   13                  14            15         16      17
            std::vector<std::string> cells; cells.reserve(20);
            std::string cur; for (char c: line) { if (c==',') { cells.push_back(cur); cur.clear(); } else { cur.push_back(c); } } cells.push_back(cur);
            auto rd = [&](int idx){ if (idx >=0 && idx < (int)cells.size()) { try { return std::stod(cells[idx]); } catch(...) { return 0.0; } } return 0.0; };
            double Pel = rd(5), Ppump = rd(17), Qspace = rd(12), Qdhw = rd(13);
            sumPel_kWh += Pel * dt_h;
            sumPump_kWh += Ppump * dt_h;
            sumQspace_kWh += Qspace * dt_h;
            sumQdhw_kWh += Qdhw * dt_h;
        }
        double elec_kWh = sumPel_kWh + sumPump_kWh;
        double heat_kWh = sumQspace_kWh + sumQdhw_kWh;
        double cost_elec = elec_kWh * std::max(0.0, cfg_.econ.elec_price_per_kWh);
        double hoursSim = (double)cfg_.time.totalSteps * cfg_.time.timeStep_s / 3600.0;
        double scale = (hoursSim > 0.0) ? (8760.0 / hoursSim) : 0.0;
        double elec_kWh_annual = elec_kWh * scale;
        double heat_kWh_annual = heat_kWh * scale;
        double r = std::max(0.0, cfg_.econ.discount_rate);
        double N = std::max(1.0, cfg_.econ.lifetime_years);
        double CRF = (r > 0.0) ? (r * std::pow(1.0 + r, N)) / (std::pow(1.0 + r, N) - 1.0) : (1.0/N);
        double LCOH_annual = (heat_kWh_annual > 1e-9) ? ((cfg_.econ.capex * CRF) + elec_kWh_annual * cfg_.econ.elec_price_per_kWh) / heat_kWh_annual : 0.0;
        std::ofstream ofs("summary.csv", std::ios::out | std::ios::trunc);
        if (ofs) {
            ofs << "elec_kWh,comp_kWh,pump_kWh,heat_kWh,space_kWh,dhw_kWh,cost_elec,capex,lifetime_y,discount,CRF,LCOH_annual\n";
            ofs << (elec_kWh) << ',' << (sumPel_kWh) << ',' << (sumPump_kWh) << ',' << (heat_kWh) << ','
                << (sumQspace_kWh) << ',' << (sumQdhw_kWh) << ',' << (cost_elec) << ','
                << (cfg_.econ.capex) << ',' << (cfg_.econ.lifetime_years) << ',' << (cfg_.econ.discount_rate) << ',' << (CRF) << ',' << (LCOH_annual) << "\n";
        }
    } catch(...) {
        // ignore errors
    }
    std::cout << std::endl;
    std::cout << "Simulation finished. Results at: " << csv_path << "\n";
    return true;
}
