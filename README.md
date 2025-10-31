# Enhance / 地埋换热与热泵系统仿真（C++ / Visual Studio）

[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)]() [![Visual Studio](https://img.shields.io/badge/IDE-Visual%20Studio%202019/2022-purple)]() [![Platform](https://img.shields.io/badge/Platform-Windows-x64)]()

地埋同轴井换热 + 热泵系统的数值仿真程序，支持半年供热/半年停运、按天气驱动的建筑负荷、可变时间步长、热泵-水箱联动、双侧流量统计与压降/泵功估算，以及批量扫参脚本（PowerShell）。

---

## 功能特性
- 季节运行：供暖季（默认 10/15～4/15），非供暖季停热，地温可自然恢复
- 变步长：前 2000 步 10 分钟；随后每步 +6 秒直至 1 小时；之后固定 2 小时
- 负荷模型：按天气温度与室内设定温度计算；当 `Tout > 26°C` 时负荷为 0（可关闭）
- 热泵 / 水箱：设定出水温 + 死区启停；去除 COP 上限（仅保留下限与防冻保护）
- 物理耦合：
  - 内/外侧对流换热（Re/Pr/Nu 动态相关式）
  - 管壁（内/外）、灌浆、土壤串联热阻
  - 顶段等效保温（长度/导热可配）、底段增强（长度可配）
- 水力/泵：摩阻 + 局部损失估算，输出压降 dP 与泵功率
- 输出：每小时一行 `results.csv`；调试 `debug.csv`；运行结束 `summary.csv`
- 扫参：`sweep.ps1` 批量运行，生成 `runs_时间戳/` 与 `summary_runs.csv`（含与实测对比）

---

## 目录结构（关键）
```
enhance/
├─ enhance.sln                   # Visual Studio 解决方案
├─ sweep.ps1                     # PowerShell 扫参脚本（在此目录运行）
├─ README.md                     # 本说明
├─ enhance/                      # 源码（VS 工程）
│  ├─ *.cpp *.h                  # SimulationController / HeatPump / Ground / Logger 等
│  ├─ weather_gansu.csv          # 甘肃天气样例
│  └─ eigen-3.4.0/               # 依赖（若存在）
└─ x64/Release/enhance.exe       # 构建产物（执行文件）
```

---

## 快速开始（Visual Studio）
1. 打开 `enhance.sln`，选择配置 `x64 | Release`
2. 右键解决方案或项目 → “重新生成”
3. 可执行文件在：`x64\Release\enhance.exe`

线程数：默认 12 线程。可用环境变量覆盖（任选其一）：
```
$env:OMP_NUM_THREADS='12'
# 或
$env:OMP_THREADS='12'
# 或按百分比
$env:OMP_THREADS_PCT='50'   # 约等于总逻辑核的 50%
```

---

## 单次运行：常用环境变量（PowerShell 示例）
```
# 天气与负荷
$env:LOAD_WEATHER_CSV='enhance\weather_gansu.csv'
$env:LOAD_INDOOR_T='26'
$env:LOAD_UA_KW_PER_K='22'
$env:LOAD_BASE_KW='0'
$env:LOAD_CUTOFF_ENABLE='1'
$env:LOAD_HEAT_CUTOFF_C='26'

# 供暖季（10/15～4/15）
$env:HEAT_SEASON_ENABLE='1'
$env:HEAT_START_MM='10'; $env:HEAT_START_DD='15'
$env:HEAT_END_MM='4';    $env:HEAT_END_DD='15'

# 热泵/水箱
$env:TANK_VOL_M3='8'
$env:TANK_SET_C='40'
$env:TANK_DB_K='2'
$env:HP_MAX_Q_OUT_KW='350'
$env:HP_EFF_CARNOT='0.65'
$env:HP_EVAP_APPROACH_K='3'
$env:HP_COND_APPROACH_K='3'
$env:HP_MIN_SRC_RETURN_C='10'
$env:HP_MOD_RANGE_C='1'
$env:HP_MAX_SRC_DT_PER_H='0'   # 0=禁用源侧降温上限

# 几何/材料/增强/保温
$env:D_OUTER_M='0.2'; $env:D_INNER_M='0.1'
$env:BORE_D_M='0.22'; $env:PIPE_THICK_M='0.008'
$env:PIPE_K_INNER='0.4'   # PE-RT II
$env:PIPE_K_OUTER='4.0'   # J55
$env:INSUL_TOP_LEN_M='100'; $env:INSUL_K_INNER='0.1'
$env:EHEP_BOTTOM_LEN_M='100'
$env:SOIL_K='2.5'; $env:SOIL_RHO='2600'; $env:SOIL_CP='1600'; $env:GROUT_K='4'

# 流量（双侧）
$env:FLOW_SRC_KGPS='10'     # 源侧 kg/s（用于地源与水力计算）
$env:FLOW_LOAD_KGPS='22.2'  # 负荷侧 kg/s（统计用）

# 运行
./x64/Release/enhance.exe
```

---

## 扫参脚本（sweep.ps1）
位置：仓库根目录。
```
powershell -NoProfile -ExecutionPolicy Bypass -File .\sweep.ps1
```
脚本会：
- 自动定位 `x64\Release\enhance.exe`
- 设置基线环境，并按脚本顶部的 `*List`（如 `$hpList`、`$uaList`、`$pipeKInnerList`）扫参
- 将每组输出归档到 `runs_时间戳/`，并生成 `summary_runs.csv`

summary_runs.csv 字段：
- `hp_max_kW, UA_kWperK, flow_src_kgps, flow_load_kgps, set_load_out_C, tank_vol_m3,
   Q_load_kWh, Q_src_kWh, P_el_kWh, P_pump_kWh, dP_kPa_avg,
   avg_Q_load_kW, avg_Q_src_kW, COP_annual, COP_measured, COP_delta,
   HP_on_hours, HP_on_measured, HP_on_delta`

> 默认与实测对比：`COP_measured=4.3`，`HP_on_measured=2890`。

---

## 输出文件说明
- `results.csv`（每小时一行）：
  - `step,date,time,T_source_out_C,T_return_C,COP,Q_out_kW,P_el_kW,Q_geo_kW,model,
     T_tank_C,HP_on,Q_space_req_kW,Q_dhw_req_kW,Q_space_served_kW,Q_dhw_served_kW,Q_unmet_kW,
     flow_src_kgps,flow_load_kgps,dP_kPa,P_pump_kW`
- `debug.csv`：热泵循环调试信息（CoolProp 调用、焓/压等）
- `summary.csv`：单次运行的能耗/热量/峰值/经济汇总
- `runs_*/results_*.csv` / `debug_*.csv`：扫参归档结果
- `runs_*/summary_runs.csv`：扫参汇总（含与实测对比项）

---

## 小贴士（性能与显示）
- 进度条按整数百分比更新并尽量单行刷新；若在 VS 的“输出”窗口出现多行堆叠，属于宿主处理回车方式不同。建议在系统 PowerShell/命令行运行以获得单行闪烁效果。
- 默认 12 线程；可用 `OMP_NUM_THREADS` / `OMP_THREADS` / `OMP_THREADS_PCT` 覆盖。

---

## Git / GitHub 提交建议
- 建议忽略：`.vs/`、`x64/`、`Debug/`、`Release/`、`*.pdb`、`*.obj`、`enhance/results.csv`、`enhance/debug.csv`、`runs_*/`、`opt_runs*/`、`enhance/summary.csv`
- VS 图形界面流程：先同步（Pull/Fetch）→ 解决冲突并提交（Commit）→ 再同步（Push）

---

## 许可
仅用于研究与教学。如需商用或二次分发，请先联系作者。