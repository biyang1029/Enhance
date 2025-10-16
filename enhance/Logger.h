#pragma once
#include <fstream>
#include <string>


class Logger {
public:
	Logger() = default;
	~Logger() { close(); }


	bool open(const std::string& path);
	void writeHour(int hour,
		double T_source_out_C,
		double T_return_C,
		double COP,
		double Q_out_kW,
		double P_el_kW,
		double Q_geo_kW);
	void close();


private:
	std::ofstream ofs_;
};