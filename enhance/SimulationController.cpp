// ========================= SimulationController.cpp =========================
#include "SimulationController.h"
#include <iostream>
#include <string>

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

    if (!logger_.open(csv_path)) {
        std::cerr << "[Error] Failed to open log file: " << csv_path << "\n";
        return false;
    }
    // 同时打开 debug 日志
    logger_.openDebug("debug.csv");

    // 回水初值
    double T_return = cfg_.fluid.inlet_T_C;
    const int    maxIter = 300;
    const double tolT    = 0.005;
    const double relax   = 0.5;

    for (int h = 0; h < cfg_.time.totalSteps; ++h) {
        // 每小时（如启用）计算天气负荷并注入
        if (cfg_.load.enable_weather) {
            static bool loaded = false;
            if (!loaded) { load_.loadCSV(cfg_.load.weather_csv, cfg_.load.column_name); loaded = true; }
            hp_.setDemand(load_.demandAt(h, cfg_.load));
        }

        double COP = 0, Q_out_kW = 0, P_el_kW = 0, Q_geo_kW = 0;
        double T_iter = T_return;
        double T_source_out = 0.0, T_next = T_iter;

        // 固定点预迭代（不推进土壤）
        for (int it = 0; it < maxIter - 1; ++it) {
            T_source_out = ground_.step(T_iter, Q_geo_kW, /*advanceSoil=*/false);
            T_next = hp_.step(T_source_out, COP, Q_out_kW, P_el_kW);
            if (std::fabs(T_next - T_iter) < tolT) break;
            T_iter = relax * T_next + (1.0 - relax) * T_iter;
        }

        // 最后一次：推进土壤并完成一步
        T_source_out = ground_.step(T_iter, Q_geo_kW, /*advanceSoil=*/true);
        double T_return_next = hp_.step(T_source_out, COP, Q_out_kW, P_el_kW);

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
            logger_.writeHour(h, T_source_out, T_return_next, COP, Q_out_kW, P_el_kW, Q_geo_kW, model);
            logger_.writeDebugHour(h, dbg.fluid, dbg.used_coolprop,
                                   dbg.T_evap_sat_K, dbg.T_cond_sat_K,
                                   dbg.P_evap_kPa, dbg.P_cond_kPa,
                                   dbg.h1, dbg.h2s, dbg.h2, dbg.h3);
        }

        T_return = T_return_next;
        printProgress(h + 1, cfg_.time.totalSteps);
    }

    logger_.close();
    std::cout << std::endl;
    std::cout << "Simulation finished. Results at: " << csv_path << "\n";
    return true;
}
