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


private:
	DataConfig cfg_;
};