#pragma once
#include <string>
#include <vector>

class RefpropWrapper {
public:
    bool load(const std::string& dll_dir);
    void setFluid(const std::string& fluid);

    // Saturation pressure [kPa] at T [K] (quality flag: 0=liq,1=vap)
    bool psat_T(double T_K, double& P_kPa, int Qflag = 1) const;

    // Enthalpy/entropy at (T [K], P [kPa])
    bool hs_from_TP(double T_K, double P_kPa, double& h_Jpkg, double& s_JpkgK) const;

    // Enthalpy and temperature from (P [kPa], s [J/kg-K])
    bool hT_from_PS(double P_kPa, double s_JpkgK, double& h_Jpkg, double& T_K) const;

private:
    void* hlib_ = nullptr; // HMODULE
    typedef void (*REFPROPdll_t)(char*, char*, char*, long&, long&, long&, double&, double&, long&, double*, double*, char*, long, long, long, long);
    REFPROPdll_t REFPROPdll_ = nullptr;
    std::string fluid_ = "R134a";
};
