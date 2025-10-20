// ========================= Logger.cpp =========================
#include "Logger.h"
#include <iomanip>


bool Logger::open(const std::string& path) {
    // Locally close any existing streams without calling close()
    if (this->ofs_.is_open()) this->ofs_.close();
    if (this->ofs_dbg_.is_open()) this->ofs_dbg_.close();
    this->ofs_.open(path, std::ios::out | std::ios::trunc);
    if (!this->ofs_) return false;
    this->ofs_ << "step,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW,model,T_tank_C,HP_on,Q_space_req_kW,Q_dhw_req_kW,Q_space_served_kW,Q_dhw_served_kW,Q_unmet_kW,flow_kgps,dP_kPa,P_pump_kW\n";
    return true;
}

bool Logger::openDebug(const std::string& path) {
    if (this->ofs_dbg_.is_open()) this->ofs_dbg_.close();
    this->ofs_dbg_.open(path, std::ios::out | std::ios::trunc);
    if (!this->ofs_dbg_) return false;
    this->ofs_dbg_ << "hour,fluid,used_coolprop,T_evap_sat_K,T_cond_sat_K,P_evap_kPa,P_cond_kPa,h1,h2s,h2,h3\n";
    return true;
}


void Logger::writeHour(int hour,
	double T_source_out_C,
	double T_return_C,
	double COP,
	double Q_out_kW,
	double P_el_kW,
	double Q_geo_kW,
	const std::string& model,
	double T_tank_C,
	int    HP_on,
	double Q_space_req_kW,
	double Q_dhw_req_kW,
	double Q_space_served_kW,
	double Q_dhw_served_kW,
	double Q_unmet_kW,
	double flow_kgps,
	double dP_kPa,
	double P_pump_kW) {
	if (!this->ofs_) return;
	this->ofs_ << hour << ','
		<< std::fixed << std::setprecision(4)
		<< T_source_out_C << ','
		<< T_return_C << ','
		<< COP << ','
		<< Q_out_kW << ','
		<< P_el_kW << ','
		<< Q_geo_kW << ',' << model << ','
		<< T_tank_C << ','
		<< HP_on << ','
		<< Q_space_req_kW << ','
		<< Q_dhw_req_kW << ','
		<< Q_space_served_kW << ','
		<< Q_dhw_served_kW << ','
		<< Q_unmet_kW << ','
		<< flow_kgps << ','
		<< dP_kPa << ','
		<< P_pump_kW << "\n";
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
    if (!this->ofs_dbg_) return;
    this->ofs_dbg_ << hour << ',' << fluid << ',' << (used_coolprop?1:0) << ','
        << std::fixed << std::setprecision(6)
        << T_evap_sat_K << ',' << T_cond_sat_K << ','
        << P_evap_kPa << ',' << P_cond_kPa << ','
        << h1 << ',' << h2s << ',' << h2 << ',' << h3 << "\n";
}


void Logger::close() {
    if (this->ofs_.is_open()) this->ofs_.close();
    if (this->ofs_dbg_.is_open()) this->ofs_dbg_.close();
}

