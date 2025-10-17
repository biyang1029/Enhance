// GroundModule.cpp
#include "GroundModule.h"
#include <algorithm>
#include <cmath>
#include <iostream>


// ---- 初始化 ----
void GroundModule::initialize(const DataConfig& cfg) {
    cfg_ = cfg;
    N_ = cfg_.well.segments;
    dz_ = cfg_.well.depth_m / static_cast<double>(N_);

    zc_.resize(N_);
    T_outer_.assign(N_, cfg_.fluid.inlet_T_C);
    T_inner_.assign(N_, cfg_.fluid.inlet_T_C);

    // 段中心深度
    for (int i = 0; i < N_; ++i) {
        zc_[i] = (i + 0.5) * dz_;
    }
    // 井壁半径（用外管外径的一半）
    double D_bore = (cfg_.well.borehole_D_m > 0.0 ? cfg_.well.borehole_D_m : cfg_.well.D_outer_m);
    r_bore_ = std::max(1e-6, 0.5 * D_bore);

    // Nz_ 与 z_（纵向与管段一致）
    Nz_ = N_;
    z_.resize(Nz_);
    for (int i = 0; i < Nz_; ++i) z_[i] = zc_[i];

    // 构建径向稳健型网格并初始化土壤温度场
    build_soil_grid_();
    init_soil_field_();

    // 初始把管内水温设为贴近壁面土壤温，减少第0小时尖峰
    double Twall0 = Tsoil_.size() ? Tsoil_[0][0] : (cfg_.T_surface_C + cfg_.geograd_C_per_m * zc_[0]);
    std::fill(T_outer_.begin(), T_outer_.end(), Twall0);
    std::fill(T_inner_.begin(), T_inner_.end(), Twall0);

    // 预分配井壁线热流
    q_line_wall_.assign(Nz_, 0.0);

}

// ---- 每小时一步：输入回水温，返回出水温；同时返回该小时取热量（kW） ----
double GroundModule::step(double T_return_C, double& Q_extracted_kW, bool advanceSoil) {


    // 井壁处土壤温度作为本小时壁温（来自上一步更新的 Tsoil_）
    std::vector<double> T_soil_local(N_);
    for (int i = 0; i < N_; ++i) {
        T_soil_local[i] = (Nr_ > 0) ? Tsoil_[0][i]
            : (cfg_.T_surface_C + cfg_.geograd_C_per_m * zc_[i]);
    }

    // 1) 外管：由井口向下（回水进入外管顶部）
    update_outer_downstream_(T_return_C, T_soil_local);

    // 2) 内管：由井底向上（内外管之间换热，先不重复计土壤）
    update_inner_upstream_(T_outer_.back(), T_outer_);

    // 3) 统一用 R′（单位长度总热阻）构造 q_line_wall_[i] (W/m)
    //    —— 这些公共量只计算一次，后面复用 —— 
    const double D_out = std::max(1e-6, cfg_.well.D_outer_m);
    const double r_po = 0.5 * D_out;
    const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
    const double r_pi = std::max(1e-6, r_po - t_p);
    const double r_bo = std::max(r_bore_, r_po + 1e-6);
    const double P_f = PI * D_out;

    const double Dh_outer = std::max(1e-6, cfg_.well.D_outer_m - cfg_.well.D_inner_m);
    const double Re_outer = cfg_.calcRe(cfg_.fluid.massFlow_kgps, Dh_outer);
    const double Pr = cfg_.calcPr();

    for (int i = 0; i < N_; ++i) {
        const double z = zc_[i];
        const double Nu = cfg_.NuFunc ? cfg_.NuFunc(Re_outer, Pr, z) : 50.0;
        const double h = Nu * cfg_.fluid.k / D_out; // W/m^2-K（外环流体侧）

        // 单位长度总热阻 R′ (K/W per meter)
        const double Rf = 1.0 / std::max(1e-9, h * P_f);
        const double Rp = std::log(r_po / r_pi) / (2.0 * PI * std::max(1e-9, cfg_.pipe.k));
        const double Rg = std::log(r_bo / r_po) / (2.0 * PI * std::max(1e-9, cfg_.grout.k));
        const double Rtot = Rf + Rp + Rg;

        const double T_wall = Tsoil_[0][i]; // 井壁处土温
        q_line_wall_[i] = (T_outer_[i] - T_wall) / Rtot; // W/m
    }

    // 4) 是否推进土壤（Neumann 边界用 q_line_wall_）
    if (advanceSoil) {
        soil_step_ADI_(cfg_.time.timeStep_s);
    }

    // 5) 用线热流积分得到 Q_extracted_kW
    double Qtot_W = 0.0;
    for (int i = 0; i < N_; ++i) Qtot_W += q_line_wall_[i] * dz_;
    Q_extracted_kW = Qtot_W / 1000.0; // kW

    // 6) 返回井口出水温（内管上行出口）
    return T_inner_.front();


}

// ---- 外管由上到下：分段推进（与土壤等效换热） ----
void GroundModule::update_outer_downstream_(double T_inlet_C, const std::vector<double>& T_soil_local) {
    const double m_dot = std::max(1e-9, cfg_.fluid.massFlow_kgps);
    const double cp = cfg_.fluid.cp;

    const double Dh_outer = std::max(1e-6, cfg_.well.D_outer_m - cfg_.well.D_inner_m);
    const double Re_outer = cfg_.calcRe(m_dot, Dh_outer);
    const double Pr = cfg_.calcPr();

    const double D_out = std::max(1e-6, cfg_.well.D_outer_m);
    const double r_po = 0.5 * D_out;
    const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
    const double r_pi = std::max(1e-6, r_po - t_p);
    const double r_bo = std::max(r_bore_, r_po + 1e-6);
    const double P_f = PI * D_out; // 流体侧换热周长

    double T_prev = T_inlet_C;

    for (int i = 0; i < N_; ++i) {
        const double z = zc_[i];
        const double Nu = cfg_.NuFunc ? cfg_.NuFunc(Re_outer, Pr, z) : (0.023 * std::pow(Re_outer, 0.8) * std::pow(Pr, 0.3));;
        const double h = Nu * cfg_.fluid.k / D_out;              // W/m2-K

        // 单位长度总热阻 R′ (K/W per meter)
        const double Rf = 1.0 / std::max(1e-9, h * P_f);
        const double Rp = std::log(r_po / r_pi) / (2.0 * PI * std::max(1e-9, cfg_.pipe.k));
        const double Rg = std::log(r_bo / r_po) / (2.0 * PI * std::max(1e-9, cfg_.grout.k));
        const double Rtot = Rf + Rp + Rg;

        const double T_wall = T_soil_local[i];
        const double q_line = (T_prev - T_wall) / Rtot;           // W/m （线热流）
        const double Q_segW = q_line * dz_;                       // W
        const double dT = -Q_segW / (m_dot * cp);

        T_outer_[i] = T_prev + dT;
        T_prev = T_outer_[i];
    }

    
}



// ---- 内管由下到上：分段推进（与外管水体换热） ----
void GroundModule::update_inner_upstream_(double T_bottom_in_C, const std::vector<double>& T_outer_local) {
    const double m_dot = std::max(1e-9, cfg_.fluid.massFlow_kgps); // kg/s
    const double cp = cfg_.fluid.cp;                             // J/kg-K
    const double D_inner = std::max(1e-6, cfg_.well.D_inner_m);
    const double P_inner = PI * D_inner;                              // 内管周长（每米长度）

    double T_prev = T_bottom_in_C; // 井底进入内管的温度（来自外管末端）

    // 由底向顶：i = N_-1 → 0
    for (int i = N_ - 1; i >= 0; --i) {
        const double T_outer = T_outer_local[i];
        const double q_io = h_inner_outer_ * (T_prev - T_outer); // W/m^2，正值：内管向外管放热
        const double Q_seg_W = q_io * (P_inner * dz_);
        const double dT = -Q_seg_W / (m_dot * cp);

        T_inner_[i] = T_prev + dT;
        T_prev = T_inner_[i];
    }

}



// ---- 估算该小时总取热量（kW）：对土壤侧热流积分 ----
double GroundModule::integrate_Q_extracted_kW_(const std::vector<double>& q_wall_Wm2) const {
    const double D_ref = std::max(1e-6, cfg_.well.D_outer_m);
    const double P_outer = PI * D_ref; // 周长（每米长度）

    double Q_total_W = 0.0;
    for (int i = 0; i < N_; ++i) {
        Q_total_W += q_wall_Wm2[i] * (P_outer * dz_); // W = W/m^2 * m^2
    }
    return Q_total_W / 1000.0; // kW
}
void GroundModule::build_soil_grid_() {
    r_.clear(); dr_.clear();
    const double rMax = cfg_.rMax_m;
    double r_acc = r_bore_;
    double dr = cfg_.dr0_m;
    while (r_acc < rMax) {
        double r_cell_center = r_acc + 0.5 * dr;
        r_.push_back(r_cell_center);
        dr_.push_back(dr);
        r_acc += dr;
        dr *= cfg_.dr_growth;
    }
    Nr_ = static_cast<int>(r_.size());
    if (Nr_ < 3) { // 防止过少
        // 回退到至少3层
        while (Nr_ < 3) {
            double back = dr_.empty() ? cfg_.dr0_m : dr_.back();
            r_.push_back((r_acc + 0.5 * back));
            dr_.push_back(back);
            r_acc += back;
            ++Nr_;
        }
    }
}

void GroundModule::init_soil_field_() {
    Tsoil_.assign(Nr_, std::vector<double>(Nz_, 0.0));
    Ttmp_.assign(Nr_, std::vector<double>(Nz_, 0.0));
    for (int iz = 0; iz < Nz_; ++iz) {
        double Tz = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[iz];
        for (int ir = 0; ir < Nr_; ++ir) {
            Tsoil_[ir][iz] = Tz;
        }
    }
}
// ADI 主入口：先 r-扫，再 z-扫（两步都是半步 dt/2 隐式）
// 物性用 cfg_.soil（你也可按半径分层扩展）
void GroundModule::soil_step_ADI_(double dt_s) {
    // 半步
    double half = 0.5 * dt_s;

    // 第一步：沿 r 扫描（每个 iz 独立解一条三对角）
    sweep_r_(half);

    // 第二步：沿 z 扫描（每个 ir 独立解一条三对角）
    sweep_z_(half);
}

// ---- r-sweep：固定 z，沿 r 解 ----
void GroundModule::sweep_r_(double dt_s) {
    const double rho = cfg_.soil.rho;
    const double cp  = cfg_.soil.cp;
    const double k   = cfg_.soil.k;

    #pragma omp parallel for schedule(static)
    for (int iz = 0; iz < Nz_; ++iz) {
        std::vector<double> a(Nr_, 0.0), b(Nr_, 0.0), c(Nr_, 0.0), d(Nr_, 0.0), x(Nr_, 0.0);

        for (int ir = 0; ir < Nr_; ++ir) Ttmp_[ir][iz] = Tsoil_[ir][iz];

        for (int ir = 0; ir < Nr_; ++ir) {
            double r   = r_[ir];
            double dr  = dr_[ir];
            double h_z = (iz == 0 || iz == Nz_ - 1)
                       ? ((Nz_ > 1) ? (z_[std::min(iz + 1, Nz_ - 1)] - z_[std::max(iz - 1, 0)]) / 2.0 : 10.0)
                       : (z_[iz + 1] - z_[iz - 1]) * 0.5;
            double area_r = 2.0 * PI * r * h_z;
            double vol    = area_r * dr;
            double cap    = rho * cp * vol / dt_s;

            double km = k, kp = k;
            double drm = (ir == 0)        ? dr : 0.5 * (dr_[ir - 1] + dr_[ir]);
            double drp = (ir == Nr_ - 1)  ? dr : 0.5 * (dr_[ir] + dr_[ir + 1]);

            double aW = (ir > 0        ? km * area_r / drm : 0.0);
            double aE = (ir < Nr_ - 1  ? kp * area_r / drp : 0.0);

            a[ir] = -aW;
            c[ir] = -aE;
            b[ir] = cap + aW + aE;
            d[ir] = cap * Tsoil_[ir][iz];
        }

        double q_line = q_line_wall_[iz];
        double q_face = q_line / (2.0 * PI * std::max(1e-6, r_bore_));
        double h_z    = (iz == 0 || iz == Nz_ - 1)
                      ? ((Nz_ > 1) ? (z_[std::min(iz + 1, Nz_ - 1)] - z_[std::max(iz - 1, 0)]) / 2.0 : 10.0)
                      : (z_[iz + 1] - z_[iz - 1]) * 0.5;
        double area_face = 2.0 * PI * std::max(1e-6, r_bore_) * h_z;
        d[0] += q_face * area_face;

        double T_far = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[iz];
        a[Nr_ - 1] = 0.0;
        c[Nr_ - 1] = 0.0;
        b[Nr_ - 1] = 1.0;
        d[Nr_ - 1] = T_far;

        thomas_solve_(a, b, c, d, x);
        for (int ir = 0; ir < Nr_; ++ir) Tsoil_[ir][iz] = x[ir];
    }
}

// ---- z-sweep：固定 r，沿 z 解 ----
void GroundModule::sweep_z_(double dt_s) {
    const double rho = cfg_.soil.rho;
    const double cp  = cfg_.soil.cp;
    const double k   = cfg_.soil.k;

    const double dz = dz_;

    #pragma omp parallel for schedule(static)
    for (int ir = 0; ir < Nr_; ++ir) {
        double r   = r_[ir];
        double dr  = dr_[ir];
        double area_z = 2.0 * PI * r * dr;
        double vol    = area_z * dz;

        std::vector<double> a(Nz_, 0.0), b(Nz_, 0.0), c(Nz_, 0.0), d(Nz_, 0.0), x(Nz_, 0.0);

        for (int iz = 0; iz < Nz_; ++iz) {
            double cap = rho * cp * vol / dt_s;

            double aS = (iz > 0         ? k * area_z / dz : 0.0);
            double aN = (iz < Nz_ - 1   ? k * area_z / dz : 0.0);

            a[iz] = -aS;
            c[iz] = -aN;
            b[iz] = cap + aS + aN;
            d[iz] = cap * Tsoil_[ir][iz];
        }

        // Dirichlet boundary at top/bottom using geothermal gradient profile
        // T(z) = T_surface + geograd * z
        if (Nz_ >= 1) {
            double T_top = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[0];
            a[0] = 0.0; c[0] = 0.0; b[0] = 1.0; d[0] = T_top;
        }
        if (Nz_ >= 2) {
            int izb = Nz_ - 1;
            double T_bot = cfg_.T_surface_C + cfg_.geograd_C_per_m * z_[izb];
            a[izb] = 0.0; c[izb] = 0.0; b[izb] = 1.0; d[izb] = T_bot;
        }

        thomas_solve_(a, b, c, d, x);
        for (int iz = 0; iz < Nz_; ++iz) Tsoil_[ir][iz] = x[iz];
    }
}

// ---- Thomas 三对角求解 ----
void GroundModule::thomas_solve_(const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    std::vector<double>& d,
    std::vector<double>& x) {
    const int n = (int)b.size();
    std::vector<double> c2(n, 0.0), d2(n, 0.0);
    c2[0] = c[0] / std::max(1e-12, b[0]);
    d2[0] = d[0] / std::max(1e-12, b[0]);
    for (int i = 1; i < n; ++i) {
        double denom = std::max(1e-12, b[i] - a[i] * c2[i - 1]);
        c2[i] = (i < n - 1) ? c[i] / denom : 0.0;
        d2[i] = (d[i] - a[i] * d2[i - 1]) / denom;
    }
    x[n - 1] = d2[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = d2[i] - c2[i] * x[i + 1];
}
