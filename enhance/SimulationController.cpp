// ========================= SimulationController.cpp =========================
#include "SimulationController.h"
#include <iostream>
// 极简进度条函数
void printProgress(int done, int total, int barWidth = 40) {
	if (done < 0) done = 0;
	if (done > total) done = total;
	double frac = static_cast<double>(done) / static_cast<double>(total);
	int filled = static_cast<int>(frac * barWidth + 0.5);

	std::cout << "\r[";
	for (int i = 0; i < barWidth; ++i) std::cout << (i < filled ? '=' : ' ');
	std::cout << "] " << int(frac * 100.0) << "% (" << done << "/" << total << ")"
		<< std::flush;
}



bool SimulationController::run(const std::string& csv_path) {
	// 初始化模块
	ground_.initialize(cfg_);
	hp_.initialize(cfg_);


	if (!logger_.open(csv_path)) {
		std::cerr << "[Error] Failed to open log file: " << csv_path << "\n";
		return false;
	}



	// 初始回水温：用配置中的外管入口温度作为第一小时回水初值
	double T_return = cfg_.fluid.inlet_T_C;
	const int    maxIter = 300;     // 每小时最多迭代 3 次
	const double tolT = 0.005;  // K，回水收敛阈值
	const double relax = 0.5;   // 欠松弛

	for (int h = 0; h < cfg_.time.totalSteps; ++h) {
		double COP = 0, Q_out_kW = 0, P_el_kW = 0, Q_geo_kW = 0;
		double T_iter = T_return;
		double T_source_out = 0.0, T_next = T_iter;

		// 预测迭代：不推进土壤（advanceSoil=false）
		for (int it = 0; it < maxIter - 1; ++it) {
			T_source_out = ground_.step(T_iter, Q_geo_kW, /*advanceSoil=*/false);
			T_next = hp_.step(T_source_out, COP, Q_out_kW, P_el_kW);

			if (std::fabs(T_next - T_iter) < tolT) break;
			T_iter = relax * T_next + (1.0 - relax) * T_iter; // 欠松弛
		}

		// 最后一遍：推进土壤（advanceSoil=true），并以收敛后的回水作为边界
		T_source_out = ground_.step(T_iter, Q_geo_kW, /*advanceSoil=*/true);
		double T_return_next = hp_.step(T_source_out, COP, Q_out_kW, P_el_kW);

		logger_.writeHour(h, T_source_out, T_return_next, COP, Q_out_kW, P_el_kW, Q_geo_kW);
		T_return = T_return_next;

		printProgress(h+1, cfg_.time.totalSteps);
	}



	logger_.close();

	std::cout << std::endl; // 结束后换行
	std::cout << "Simulation finished. Results at: " << csv_path << "\n";
	return true;
}