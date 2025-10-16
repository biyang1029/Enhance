// DataConfig.h
#pragma once
#include <functional>
#include <string>
#include <cmath>

// 圆周率常量（全局可用）
constexpr double PI = 3.14159265358979323846;

// —— 流体、材料、几何与时间配置 —— //
struct FluidProperty {
    double rho = 997.0;   // kg/m3
    double cp = 4186.0;  // J/kg-K
    double k = 0.6;     // W/m-K
    double mu = 0.001;   // Pa·s
    // 入口与流量
    double inlet_T_C = 15.0; // 井口外管入口温度
    double massFlow_kgps = 10;  // 质量流量 (kg/s)
};

struct MaterialProperty {
    double k = 2.0;     // W/m-K
    double rho = 1800.0;  // kg/m3
    double cp = 900.0;   // J/kg-K
};

struct WellGeometry {
    double depth_m = 2500.0;
    int    segments = 250;      // 纵向段数
    double D_outer_m = 0.20;     // 外管直径（示例）
    double D_inner_m = 0.10;     // 内管直径（示例）
    double borehole_D_m = 0.22;           // 井孔直径（回填外半径）
    double pipe_wall_thickness_m = 0.008; // 外管壁厚
};

struct SoilMesh {
    double rMax_m = 100.0;
    double dr0_m = 0.01;     // 第一层网格厚度
    double growth = 1.2;      // 径向增长率
};

struct TimeControl {
    double timeStep_s = 600.0;     // 1小时
    int    totalSteps = 6*24*10;         // 先跑一天
    int    maxIter = 500;         // 管内迭代上限
    double residualTol = 0.00001;       // 收敛阈值(℃)
};


// —— 热泵/负荷的最小占位（后续细化） —— //
struct HeatPumpConfig {

    double defaultCOP = 3.0;          // 可保留备用
    double Q_demand_kW = 350.0;       // 目标供热功率（会被负反馈缩放）
    double target_T_load_out_C = 45.0;// 负荷侧目标供水温（影响 COP 估计）
    double freeze_guard_C = 1.0;      // 源侧回水最低保护温度
    double mod_range_C = 10.0;        // 低温降载的线性范围（T_freeze ~ T_freeze+range）
    double max_source_dT_per_hour = 2.0; // 每小时源侧最大降温（K/h）限制
    double evap_approach_K = 3.0;   // evaporator approach (K)
    double cond_approach_K = 5.0;   // condenser approach (K)
    double eff_carnot      = 0.50;  // fraction of Carnot COP


    // Refrigerant & cycle parameters
    std::string fluid = "R134a";   // default working fluid
    double superheat_K = 3.0;       // evaporator outlet superheat
    double subcool_K   = 5.0;       // condenser outlet subcooling
    double eta_isentropic = 0.70;   // compressor isentropic efficiency
    bool   use_refprop = true;      // allow disabling REFPROP at runtime
    std::string refprop_dll_dir = "D:/yangb/soft/refpro"; // REFPROP DLL directory (empty for fallback)
};

// —— 总配置 —— //
struct DataConfig {
    // 几何/网格/时间
    WellGeometry well;
    SoilMesh     mesh;
    TimeControl  time;
    // 地温初始剖面
    double T_surface_C = 6.0;    // 地表温度 (°C)
    double geograd_C_per_m = 0.03;   // 地温梯度 (°C/m) = 3°C/100m

    // 土壤域与网格（稳健型）
    double rMax_m = 200.0;       // 土壤域半径 (m)
    double dr0_m = 0.005;       // 第一层径向厚度 (m)
    double dr_growth = 1.05;        // 径向增长率

    // 物性
    FluidProperty    fluid;
    MaterialProperty soil;     // 原状土
    MaterialProperty grout;    // 回填
    MaterialProperty pipe;     // 管壁

    // 热泵/负荷
    HeatPumpConfig   hp;

    // Nu 计算策略：Nu = f(Re, Pr, z)
    // 你可以在 cpp 里设置默认策略：z >= 1500m 使用增强段
    std::function<double(double, double, double)> NuFunc;

    // 工具函数：基于当前配置计算 Re / Pr（可在管段里用）
    inline double calcRe(double m_dot_kgps, double D_m) const {
        // Re = 4 m_dot / (pi D mu)  —— 假设密度接近常数，使用体积/剪切量纲近似
        return (4.0 * m_dot_kgps) / (PI * D_m * fluid.mu);
    }
    inline double calcPr() const {
        // Pr = cp * mu / k
        return (fluid.cp * fluid.mu) / fluid.k;
    }
};

// 全局配置（方便各模块直接 include 使用）
extern DataConfig CFG;
