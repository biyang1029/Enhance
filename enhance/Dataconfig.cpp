// DataConfig.cpp
#include "DataConfig.h"
#include <algorithm>
#include <cmath>

DataConfig CFG;  // 全局实例

// 可在此处自定义/覆盖默认参数与 Nu 策略
// 你也可以在 main 或 SimulationController 初始化时手动设置。
static struct _BootInit {
    _BootInit() {
        // 示例材料（可按需调整）
        CFG.soil = { 2.0, 1800.0, 900.0 };   // k, rho, cp
        CFG.grout = { 1.5, 2000.0, 1000.0 };
        CFG.pipe = { 15.0, 7800.0, 500.0 };

        // 默认 Nu 策略：z >= 1500m 采用“增强段”经验式，否则 Dittus-Boelter
        CFG.NuFunc = [](double Re, double Pr, double z_m) -> double {
            if (z_m >= 1500.0) {
                // 增强段（示例系数，可后续替换/查表）
                return 0.035 * std::pow(Re, 0.85) * std::pow(Pr, 0.40);
            }
            else {
                // Dittus-Boelter（加热：n=0.4；冷却：n=0.3，按需替换）
                return 0.023 * std::pow(Re, 0.80) * std::pow(Pr, 0.40);
            }
            };

        // 你也可以在这里调整时间、负荷等默认值：
        // CFG.time.totalSteps = 8760; // 跑一年
        // CFG.hp.Q_demand_kW  = 80.0;
    }
} _bootInitOnce;  // 程序加载时执行一次
