// ========================= Logger.cpp =========================
#include "Logger.h"
#include <iomanip>


bool Logger::open(const std::string& path) {
	close();
	ofs_.open(path, std::ios::out | std::ios::trunc);
	if (!ofs_) return false;
	ofs_ << "hour,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW\n";
	return true;
}


void Logger::writeHour(int hour,
	double T_source_out_C,
	double T_return_C,
	double COP,
	double Q_out_kW,
	double P_el_kW,
	double Q_geo_kW) {
	if (!ofs_) return;
	ofs_ << hour << ','
		<< std::fixed << std::setprecision(4)
		<< T_source_out_C << ','
		<< T_return_C << ','
		<< COP << ','
		<< Q_out_kW << ','
		<< P_el_kW << ','
		<< Q_geo_kW << "\n";
}


void Logger::close() {
	if (ofs_.is_open()) ofs_.close();
}