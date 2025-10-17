#pragma once
#include "DataConfig.h"


class HeatPumpModule {
public:
    HeatPumpModule() = default;
    void initialize(const DataConfig& cfg) { cfg_ = cfg; }


	// 输入：地源侧出水温（井口内管出口温度，供给热泵蒸发器）
	// 输出：回到地热回路的回水温（进入外管顶部的温度）
	// 通过引用返回：COP、供热量Q_out(kW)、电功率P_el(kW)
    double step(double T_source_out_C, double& COP, double& Q_out_kW, double& P_el_kW);

    struct HPDebug {
        bool used_coolprop = false;
        std::string fluid;
        double P_evap_kPa = 0.0, P_cond_kPa = 0.0;
        double T_evap_sat_K = 0.0, T_cond_sat_K = 0.0;
        double h1 = 0.0, h2s = 0.0, h2 = 0.0, h3 = 0.0;
    };
    void setDemand(double q_kW) { cfg_.hp.Q_demand_kW = q_kW; }
    const HPDebug& lastDebug() const { return last_debug_; }

private:
    DataConfig cfg_;
    HPDebug    last_debug_;
};
