// ========================= Logger.cpp =========================
#include "Logger.h"
#include <iomanip>


bool Logger::open(const std::string& path) {
    close();
    ofs_.open(path, std::ios::out | std::ios::trunc);
    if (!ofs_) return false;
    ofs_ << "hour,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW,model\\n";
    return true;
}

bool Logger::openDebug(const std::string& path) {
    if (ofs_dbg_.is_open()) ofs_dbg_.close();
    ofs_dbg_.open(path, std::ios::out | std::ios::trunc);
    if (!ofs_dbg_) return false;
    ofs_dbg_ << "hour,fluid,used_coolprop,T_evap_sat_K,T_cond_sat_K,P_evap_kPa,P_cond_kPa,h1,h2s,h2,h3\n";
    return true;
}


void Logger::writeHour(int hour,
	double T_source_out_C,
	double T_return_C,
	double COP,
	double Q_out_kW,
	double P_el_kW,
	double Q_geo_kW,
	const std::string& model) {
	if (!ofs_) return;
	ofs_ << hour << ','
		<< std::fixed << std::setprecision(4)
		<< T_source_out_C << ','
		<< T_return_C << ','
		<< COP << ','
		<< Q_out_kW << ','
		<< P_el_kW << ','
		<< Q_geo_kW << "," << model << "\\n";
}

void Logger::writeDebugHour(int hour,
    const std::string& fluid,
    bool used_coolprop,
    double T_evap_sat_K,
    double T_cond_sat_K,
    double P_evap_kPa,
    double P_cond_kPa,
    double h1,
    double h2s,
    double h2,
    double h3) {
    if (!ofs_dbg_) return;
    ofs_dbg_ << hour << ',' << fluid << ',' << (used_coolprop?1:0) << ','
        << std::fixed << std::setprecision(6)
        << T_evap_sat_K << ',' << T_cond_sat_K << ','
        << P_evap_kPa << ',' << P_cond_kPa << ','
        << h1 << ',' << h2s << ',' << h2 << ',' << h3 << "\n";
}


void Logger::close() {
    if (ofs_.is_open()) ofs_.close();
    if (ofs_dbg_.is_open()) ofs_dbg_.close();
}

