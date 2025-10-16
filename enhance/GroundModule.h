// GroundModule.h
#pragma once
#include "DataConfig.h"
#include <vector>

// 说明：这是“科研版最小骨架”。
// step() 内部先用“等效热阻 + 分段推进”的快速近似，先跑通闭环；
// 之后你可以把已有 Pipe/Soil 的细化实现迁入 .cpp（接口不变）。

class GroundModule {
public:
    GroundModule() = default;

    // 初始化：生成纵向段中心深度、预分配缓存等
    void initialize(const DataConfig& cfg);

    // 输入：回水温（内管出口 → 外管入口的回到地热侧的水温）
    // 输出：该小时的井口出水温（内管出口）；并通过引用返回从地层提取的热量（kW）
    double step(double T_return_C, double& Q_extracted_kW, bool advanceSoil = true);

    int    getNz() const { return Nz_; }
    double getDz() const { return dz_; }
    double getWallT(int iz) const { return (Nr_ > 0 && iz >= 0 && iz < Nz_) ? Tsoil_[0][iz] : 0.0; }

private:
    // —— 本地缓存 —— //
    DataConfig cfg_;
    int    N_ = 0;           // 段数
    double dz_ = 0.0;        // 段长
    std::vector<double> zc_; // 段中心深度（m）

    // 分段水温缓存（外管下行、内管上行）
    std::vector<double> T_outer_;  // size N_
    std::vector<double> T_inner_;  // size N_

    // 等效换热参数（先做成简化，后续替换成更精细的求解）
    // h_soil_outer：外管对土壤/回填等效换热系数（W/m2-K）
    // h_inner_outer：内外管之间等效换热（W/m2-K）
    double h_inner_outer_ = 10.0; // 先给常数；后续可基于 Nu 计算

    // —— 内部帮助函数 —— //
    // 计算某段外管的换热系数（对土壤等效），可用 Nu(Re,Pr,z) 近似

    // 计算“沿程”温度推进：T_out = T_in + (q * A * dz) / (m_dot * cp)；q = h * (T_fluid - T_wall)
    // 这里的 T_wall 先用“近似的土壤局部温度”（比如初温或上一小时更新值），先跑闭环
    // 后续替换为隐式二维土壤求解器提供的壁面温度/热流。
    void update_outer_downstream_(double T_inlet_C, const std::vector<double>& T_soil_local);
    void update_inner_upstream_(double T_bottom_in_C, const std::vector<double>& T_outer_local);

    // 估算这一小时的总取热量（kW），按段积分：sum(q_wall * 周长 * dz) / 1000
    double integrate_Q_extracted_kW_(const std::vector<double>& q_wall_Wm2) const;
    // GroundModule.h 内 class GroundModule { ... } 的 private: 增补以下成员

// —— 土壤二维网格与场（轴对称） ——
// 网格尺寸
    int    Nr_ = 0;                       // 径向层数
    int    Nz_ = 0;                       // 纵向层数（= N_）
    double r_bore_ = 0.0;                 // 井壁半径 = D_outer/2

    // 网格节点（cell-centered）
    std::vector<double> r_;               // r 方向单元中心半径 (size Nr_)
    std::vector<double> dr_;              // r 方向单元厚度 (size Nr_)
    std::vector<double> z_;               // z 方向单元中心 (size Nz_) —— 已有 zc_ 可重用

    // 温度场：Tsoil_[ir][iz]
    std::vector<std::vector<double>> Tsoil_;

    // ADI 暂存
    std::vector<std::vector<double>> Ttmp_;

    // 井壁线热流（W/m），由外管计算得到（每段一个）
    std::vector<double> q_line_wall_;     // size Nz_

    // —— 土壤网格/场 初始化与时间步推进 ——

    // 生成稳健型径向网格（dr0=0.005, growth=1.05，到 rMax=200m）
    void build_soil_grid_();

    // 用地温初始剖面填充 Tsoil：T = T_surface + geograd * z
    void init_soil_field_();

    // 每小时一次：ADI 两步推进（先 r-扫后 z-扫或反之）
    void soil_step_ADI_(double dt_s);

    // ADI 步骤 1：固定 z，沿 r 做三对角求解（Neumann 井壁热流、外边界绝热）
    void sweep_r_(double dt_s);

    // ADI 步骤 2：固定 r，沿 z 做三对角求解（顶底绝热，可改 Dirichlet）
    void sweep_z_(double dt_s);

    // Thomas 解三对角
    void thomas_solve_(const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        std::vector<double>& d,
        std::vector<double>& x);

};