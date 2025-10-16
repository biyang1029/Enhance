// DataConfig.cpp
#include "DataConfig.h"
#include <algorithm>
#include <cmath>

DataConfig CFG;  // ȫ��ʵ��

// ���ڴ˴��Զ���/����Ĭ�ϲ����� Nu ����
// ��Ҳ������ main �� SimulationController ��ʼ��ʱ�ֶ����á�
static struct _BootInit {
    _BootInit() {
        // ʾ�����ϣ��ɰ��������
        CFG.soil = { 2.0, 1800.0, 900.0 };   // k, rho, cp
        CFG.grout = { 1.5, 2000.0, 1000.0 };
        CFG.pipe = { 15.0, 7800.0, 500.0 };

        // Ĭ�� Nu ���ԣ�z >= 1500m ���á���ǿ�Ρ�����ʽ������ Dittus-Boelter
        CFG.NuFunc = [](double Re, double Pr, double z_m) -> double {
            if (z_m >= 1500.0) {
                // ��ǿ�Σ�ʾ��ϵ�����ɺ����滻/���
                return 0.035 * std::pow(Re, 0.85) * std::pow(Pr, 0.40);
            }
            else {
                // Dittus-Boelter�����ȣ�n=0.4����ȴ��n=0.3�������滻��
                return 0.023 * std::pow(Re, 0.80) * std::pow(Pr, 0.40);
            }
            };

        // ��Ҳ�������������ʱ�䡢���ɵ�Ĭ��ֵ��
        // CFG.time.totalSteps = 8760; // ��һ��
        // CFG.hp.Q_demand_kW  = 80.0;
    }
} _bootInitOnce;  // �������ʱִ��һ��
