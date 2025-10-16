// ========================= SimulationController.h =========================
#pragma once
#include <string>
#include "DataConfig.h"
#include "GroundModule.h"
#include "HeatPumpModule.h"
#include "Logger.h"


class SimulationController {
public:
	explicit SimulationController(const DataConfig& cfg) : cfg_(cfg) {}
	bool run(const std::string& csv_path = "results.csv");


private:
	DataConfig cfg_;
	GroundModule ground_;
	HeatPumpModule hp_;
	Logger logger_;
};