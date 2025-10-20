#pragma once
#include <fstream>
#include <string>


class Logger {
public:
	Logger() = default;
	~Logger() { close(); }


    bool open(const std::string& path);
    bool openDebug(const std::string& path);
    void writeHour(int hour,
        double T_source_out_C,
        double T_return_C,
        double COP,
        double Q_out_kW,
        double P_el_kW,
        double Q_geo_kW,
        const std::string& model,
        const std::string& date_ymd,
        const std::string& time_hm,
        double T_tank_C,
        int    HP_on,
        double Q_space_req_kW,
        double Q_dhw_req_kW,
        double Q_space_served_kW,
        double Q_dhw_served_kW,
        double Q_unmet_kW,
        double flow_src_kgps,
        double flow_load_kgps,
        double dP_kPa,
        double P_pump_kW);
    void writeDebugHour(int hour,
        const std::string& fluid,
        bool used_coolprop,
        double T_evap_sat_K,
        double T_cond_sat_K,
        double P_evap_kPa,
        double P_cond_kPa,
        double h1,
        double h2s,
        double h2,
        double h3);
    void close();


private:
    std::ofstream ofs_;
    std::ofstream ofs_dbg_;
};


