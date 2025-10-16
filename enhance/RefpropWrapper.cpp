#include "RefpropWrapper.h"
#include <cstring>
#include <algorithm>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#include <windows.h>
#endif

static inline void strpad(char* dst, size_t cap, const std::string& s){
    size_t n = (std::min)(cap-1, s.size());
    memcpy(dst, s.data(), n);
    dst[n] = '\0';
}

bool RefpropWrapper::load(const std::string& dll_dir){
#ifdef _WIN32
    std::string dll = dll_dir;
    if (!dll.empty() && (dll.back()!='\\' && dll.back()!='/')) dll += "\\";
    std::string path = dll + "REFPROP64.dll";
    HMODULE h = LoadLibraryA(path.c_str());
    if (!h){
        path = dll + "REFPROP.dll"; // fallback name
        h = LoadLibraryA(path.c_str());
    }
    if (!h) return false;
    hlib_ = h;
    REFPROPdll_ = (REFPROPdll_t) GetProcAddress(h, "REFPROPdll");
    return REFPROPdll_ != nullptr;
#else
    return false;
#endif
}

void RefpropWrapper::setFluid(const std::string& fluid){ fluid_ = fluid; }

bool RefpropWrapper::psat_T(double T_K, double& P_kPa, int Qflag){
    if (!REFPROPdll_) return false;
    // REFPROPdll(hFld,hIn,hOut,iUnits,iMass,iFlag,a,b,iErr,x,y,hErr,l1,l2,l3,lErr)
    char hFld[10000]={0}; strpad(hFld,sizeof(hFld),fluid_);
    char hIn[4]={0};  strpad(hIn,sizeof(hIn),"TQ");
    char hOut[8]={0}; strpad(hOut,sizeof(hOut),"P");
    long iUnits=21, iMass=1, iFlag=0, iErr=0;
    double a=T_K, b=(double)Qflag;
    double x[20]={0}; x[0]=1.0;
    double y[200]={0};
    char hErr[256]={0};
    REFPROPdll_(hFld,hIn,hOut,iUnits,iMass,iFlag,a,b,iErr,x,y,hErr,(long)strlen(hFld),(long)strlen(hIn),(long)strlen(hOut),(long)sizeof(hErr));
    if (iErr!=0) return false;
    P_kPa = y[0];
    return true;
}

bool RefpropWrapper::hs_from_TP(double T_K, double P_kPa, double& h_Jpkg, double& s_JpkgK){
    if (!REFPROPdll_) return false;
    char hFld[10000]={0}; strpad(hFld,sizeof(hFld),fluid_);
    char hIn[4]={0};  strpad(hIn,sizeof(hIn),"TP");
    char hOut[8]={0}; strpad(hOut,sizeof(hOut),"H;S");
    long iUnits=21, iMass=1, iFlag=0, iErr=0;
    double a=T_K, b=P_kPa;
    double x[20]={0}; x[0]=1.0;
    double y[200]={0};
    char hErr[256]={0};
    REFPROPdll_(hFld,hIn,hOut,iUnits,iMass,iFlag,a,b,iErr,x,y,hErr,(long)strlen(hFld),(long)strlen(hIn),(long)strlen(hOut),(long)sizeof(hErr));
    if (iErr!=0) return false;
    h_Jpkg = y[0];
    s_JpkgK = y[1];
    return true;
}

bool RefpropWrapper::hT_from_PS(double P_kPa, double s_JpkgK, double& h_Jpkg, double& T_K){
    if (!REFPROPdll_) return false;
    char hFld[10000]={0}; strpad(hFld,sizeof(hFld),fluid_);
    char hIn[4]={0};  strpad(hIn,sizeof(hIn),"PS");
    char hOut[8]={0}; strpad(hOut,sizeof(hOut),"H;T");
    long iUnits=21, iMass=1, iFlag=0, iErr=0;
    double a=P_kPa, b=s_JpkgK;
    double x[20]={0}; x[0]=1.0;
    double y[200]={0};
    char hErr[256]={0};
    REFPROPdll_(hFld,hIn,hOut,iUnits,iMass,iFlag,a,b,iErr,x,y,hErr,(long)strlen(hFld),(long)strlen(hIn),(long)strlen(hOut),(long)sizeof(hErr));
    if (iErr!=0) return false;
    h_Jpkg = y[0];
    T_K    = y[1];
    return true;
}