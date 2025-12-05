// DataConfig.cpp
#include "DataConfig.h"
#include <algorithm>
#include <cmath>

DataConfig CFG;

// Boot-time init for defaults and Nu correlation
static struct _BootInit {
    _BootInit() {
        // Default Nu model: enhanced when z>=1500m; else Dittus-Boelter
        CFG.NuFunc = [](double Re, double Pr, double z_m) -> double {
            if (z_m >= 1500.0) {
                return 0.035 * std::pow(Re, 0.85) * std::pow(Pr, 0.40);
            } else {
                return 0.023 * std::pow(Re, 0.80) * std::pow(Pr, 0.40);
            }
        };
    }
} _bootInitOnce;
