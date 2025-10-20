// GroundModule.cpp
#include "GroundModule.h"
#include <algorithm>
#include <cmath>
#include <iostream>


// ---- ????? ----
void GroundModule::initialize(const DataConfig& cfg) {
    cfg_ = cfg;
    N_ = cfg_.well.segments;
    dz_ = cfg_.well.depth_m / static_cast<double>(N_);

    zc_.resize(N_);
    T_outer_.assign(N_, cfg_.fluid.inlet_T_C);
    T_inner_.assign(N_, cfg_.fluid.inlet_T_C);

    // ?????????
    for (int i = 0; i < N_; ++i) {
        zc_[i] = (i + 0.5) * dz_;
    }
    // ???????????????????
    double D_bore = (cfg_.well.borehole_D_m > 0.0 ? cfg_.well.borehole_D_m : cfg_.well.D_outer_m);
    r_bore_ = std::max(1e-6, 0.5 * D_bore);

    // Nz_ ?? z_????????????��?
    Nz_ = N_;
    z_.resize(Nz_);
    for (int i = 0; i < Nz_; ++i) z_[i] = zc_[i];

    // ??????????????????????????????
    build_soil_grid_();
    init_soil_field_();

    // ???????????????????????????��??????0��????
    double Twall0 = Tsoil_.size() ? Tsoil_[0][0] : (cfg_.T_surface_C + cfg_.geograd_C_per_m * zc_[0]);
    std::fill(T_outer_.begin(), T_outer_.end(), Twall0);
    std::fill(T_inner_.begin(), T_inner_.end(), Twall0);

    // ????��????????
    q_line_wall_.assign(Nz_, 0.0);

}

// ---- ?��?????????????��????????��????????��????????kW?? ----
double GroundModule::step(double T_return_C, double& Q_extracted_kW, double dt_s, bool advanceSoil) {


    // ?????????????????��????��?????????????��? Tsoil_??
    std::vector<double> T_soil_local(N_);
    for (int i = 0; i < N_; ++i) {
        T_soil_local[i] = (Nr_ > 0) ? Tsoil_[0][i]
            : (cfg_.T_surface_C + cfg_.geograd_C_per_m * zc_[i]);
    }

    // 1) ????????????��????????????????
    update_outer_downstream_(T_return_C, T_soil_local);

    // 2) ?????????????????????��?????????????????
    update_inner_upstream_(T_outer_.back(), T_outer_);

    // 3) ???? R????��?????????�h???? q_line_wall_[i] (W/m)
    //    ???? ??��?????????????��????��?? ???? 
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
        const double h = Nu * cfg_.fluid.k / D_out; // W/m^2-K?????????

        // ??��?????????? R?? (K/W per meter)
        const double Rf = 1.0 / std::max(1e-9, h * P_f);
        const double k_outer = (cfg_.pipe_k_outer > 0.0 ? cfg_.pipe_k_outer : cfg_.pipe.k);
        const double Rp = std::log(r_po / r_pi) / (2.0 * PI * std::max(1e-9, k_outer));
        const double Rg = std::log(r_bo / r_po) / (2.0 * PI * std::max(1e-9, cfg_.grout.k));
        const double Rtot = Rf + Rp + Rg;

        const double T_wall = Tsoil_[0][i]; // ?????????
        q_line_wall_[i] = (T_outer_[i] - T_wall) / Rtot; // W/m
    }

    // 4) ????????????Neumann ????? q_line_wall_??
    if (advanceSoil) {
        soil_step_ADI_(dt_s);
    }

    // 5) ?????????????? Q_extracted_kW
    double Qtot_W = 0.0;
    for (int i = 0; i < N_; ++i) Qtot_W += q_line_wall_[i] * dz_;
    Q_extracted_kW = Qtot_W / 1000.0; // kW

    // 6) ??????????��???????��????
    return T_inner_.front();


}

// ---- ?????????��?????????????????��????? ----
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
    const double P_f = PI * D_out; // ??????????

    double T_prev = T_inlet_C;

    for (int i = 0; i < N_; ++i) {
        const double z = zc_[i];
        const double Nu = cfg_.NuFunc ? cfg_.NuFunc(Re_outer, Pr, z) : (0.023 * std::pow(Re_outer, 0.8) * std::pow(Pr, 0.3));;
        const double h = Nu * cfg_.fluid.k / D_out;              // W/m2-K

        // resistances per unit length
        const double Rf = 1.0 / std::max(1e-9, h * P_f);
        const double k_outer = (cfg_.pipe_k_outer > 0.0 ? cfg_.pipe_k_outer : cfg_.pipe.k);
        const double Rp = std::log(r_po / r_pi) / (2.0 * PI * std::max(1e-9, k_outer));
        const double Rg = std::log(r_bo / r_po) / (2.0 * PI * std::max(1e-9, cfg_.grout.k));
        const double Rtot = Rf + Rp + Rg;

        const double T_wall = T_soil_local[i];
        const double q_prime = (T_prev - T_wall) / std::max(1e-12, Rtot); // W/m
        const double dT = - (q_prime * dz_) / (m_dot * cp);
        T_outer_[i] = T_prev + dT;
        T_prev = T_outer_[i];
    }

    
}



// ---- ??????��???????????????????�h??? ----
void GroundModule::update_inner_upstream_(double T_bottom_in_C, const std::vector<double>& T_outer_local) {
    const double m_dot = std::max(1e-9, cfg_.fluid.massFlow_kgps); // kg/s
    const double cp    = cfg_.fluid.cp;                            // J/kg-K
    const double rho   = cfg_.fluid.rho;
    const double mu    = std::max(1e-9, cfg_.fluid.mu);
    const double kf    = std::max(1e-9, cfg_.fluid.k);             // W/m-K

    // Geometry for inner/annulus coupling
    const double t_p = std::max(1e-6, cfg_.well.pipe_wall_thickness_m);
    const double D_o_out = std::max(1e-6, cfg_.well.D_outer_m);
    const double D_o_inn = std::max(1e-6, cfg_.well.D_inner_m);
    const double D_o_inner_wall = std::max(1e-6, D_o_out - 2.0 * t_p);
    const double D_i_inner = std::max(1e-6, D_o_inn - 2.0 * t_p);
    const double r_i_inner = 0.5 * D_i_inner;
    const double r_i_outer = 0.5 * D_o_inn;

    // Hydraulics cross-sections
    const double Dh_ann = std::max(1e-6, std::fabs(D_o_inner_wall - D_o_inn));
    const double A_ann  = std::max(1e-12, 0.25 * PI * (D_o_inner_wall*D_o_inner_wall - D_o_inn*D_o_inn));
    const double A_in   = std::max(1e-12, 0.25 * PI * (D_i_inner*D_i_inner));

    auto Nu_turb = [](double Re, double Pr){ return 0.023 * std::pow(Re, 0.8) * std::pow(Pr, 0.4); }; // Dittus–Boelter
    auto Nu_lam  = [](double /*Re*/, double /*Pr*/){ return 4.36; };                                   // laminar, const-q
    auto blendNu = [&](double Re, double Pr){
        if (Re < 2000.0) return Nu_lam(Re,Pr);
        if (Re > 3000.0) return Nu_turb(Re,Pr);
        double a = (Re - 2000.0) / 1000.0; // linear blend in transition
        return (1.0 - a) * Nu_lam(Re,Pr) + a * Nu_turb(Re,Pr);
    };

    const double Pr = cfg_.calcPr();
    const double h_min = 5.0; // W/m2-K

    double T_prev = T_bottom_in_C;
    for (int i = N_ - 1; i >= 0; --i) {
        const double T_outer = T_outer_local[i];

        // Reynolds and Nusselt
        const double v_in  = m_dot / (rho * A_in);
        const double v_ann = m_dot / (rho * A_ann);
        const double Re_in  = rho * v_in  * std::max(1e-6, D_i_inner) / mu;
        const double Re_ann = rho * v_ann * std::max(1e-6, Dh_ann)    / mu;
        const double Nu_in  = blendNu(std::max(1e-6, Re_in),  std::max(1e-6, Pr));
        const double Nu_ann = blendNu(std::max(1e-6, Re_ann), std::max(1e-6, Pr));
        const double h_i = std::max(h_min, Nu_in  * kf / std::max(1e-6, D_i_inner));
        const double h_o = std::max(h_min, Nu_ann * kf / std::max(1e-6, Dh_ann));

        // Series resistances per unit length
        const double R_i = 1.0 / std::max(1e-9, h_i * 2.0 * PI * std::max(1e-9, r_i_inner));
        const double z = zc_[i];
        double kin = (cfg_.pipe_k_inner > 0.0 ? cfg_.pipe_k_inner : cfg_.pipe.k);
        if (cfg_.insul.enable && z <= std::max(0.0, cfg_.insul.top_len_m)) {
            kin = std::max(1e-4, cfg_.insul.k_inner);
        }
        const double R_w = std::log(std::max(1e-6, r_i_outer) / std::max(1e-6, r_i_inner)) / (2.0 * PI * std::max(1e-9, kin));
        const double R_o = 1.0 / std::max(1e-9, h_o * 2.0 * PI * std::max(1e-9, r_i_outer));
        const double U_per_len = 1.0 / std::max(1e-12, R_i + R_w + R_o);

        // Heat per unit length and temperature drop along dz
        const double q_prime_Wpm = U_per_len * (T_prev - T_outer);
        const double dT_raw = - (q_prime_Wpm * dz_) / (m_dot * cp);
        double dT = dT_raw; if (dT > 10.0) dT = 10.0; else if (dT < -10.0) dT = -10.0;
        T_inner_[i] = T_prev + dT;
        T_prev = T_inner_[i];
    }
}



// ---- ?????��??????????kW???????????????????? ----
double GroundModule::integrate_Q_extracted_kW_(const std::vector<double>& q_wall_Wm2) const {
    const double D_ref = std::max(1e-6, cfg_.well.D_outer_m);
    const double P_outer = PI * D_ref; // ????????????

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
    if (Nr_ < 3) { // ???????
        // ?????????3??
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
// ADI ???????? r-????? z-???????????? dt/2 ?????
void GroundModule::soil_step_ADI_(double dt_s) {
    if (Nr_ <= 0 || Nz_ <= 0) return;
    double half = std::max(1e-9, 0.5 * dt_s);

    // ????????? r ??��??? iz ???????????????
    sweep_r_(half);

    // ????????? z ??��??? ir ???????????????
    sweep_z_(half);
}

// ---- r-sweep????? z???? r ?? ----
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

// ---- z-sweep????? r???? z ?? ----
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

// ---- Thomas ???????? ----
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



