#!/usr/bin/env python3
"""
Sweep bottom enhancement length and compute multi-objective Pareto front.

Objectives (to minimize):
 - LCOH_annual (currency/kWh)
 - elec_kWh (total electricity = compressor + pump)
 - unmet_kWh (constraint-style: prefer 0; included as objective)

Derived metrics reported:
 - comp_kWh, pump_kWh, heat_kWh, SCOP_sys = heat_kWh / (comp_kWh + pump_kWh)

Usage (PowerShell):
  python optimize_enhancement.py \
    --exe "D:\\...\\enhance\\x64\\Debug\\enhance.exe" \
    --weather "D:\\...\\weather_year.csv" \
    --length-grid 0:800:9 \
    --years 1 \
    --econ-elec 0.75 --econ-capex 0 --econ-life 20 --econ-dr 0.08

Notes:
 - Enhancement is forced to bottom-only (z_start = depth - L, z_end = depth).
 - Depth is read from DataConfig.h unless --depth provided.
 - Each run writes results into a per-run folder under ./runs.
 - Requires the C++ binary to be already built.
"""

import argparse, csv, math, os, re, shutil, subprocess, sys, time
from pathlib import Path


def read_depth_from_header(root: Path) -> float:
    hdr = root / 'enhance' / 'enhance' / 'DataConfig.h'
    if not hdr.exists():
        # backtrack one level (if script in tools)
        hdr = root.parent / 'enhance' / 'enhance' / 'DataConfig.h'
    depth = 2500.0
    try:
        txt = hdr.read_text(encoding='utf-8', errors='ignore')
        m = re.search(r'double\s+depth_m\s*=\s*([0-9eE\.\+\-]+)', txt)
        if m:
            depth = float(m.group(1))
    except Exception:
        pass
    return depth


def parse_grid(spec: str):
    # format: start:end:count
    if ':' in spec:
        a, b, n = spec.split(':')
        a = float(a); b = float(b); n = int(n)
        if n < 2:
            return [a]
        step = (b - a) / (n - 1)
        return [a + i*step for i in range(n)]
    else:
        # comma list
        return [float(x) for x in spec.split(',') if x.strip()]


def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)


def run_case(exe: Path, workdir: Path, weather: Path, years: int, timestep_s: int,
             depth_m: float, L_enh_m: float,
             econ_elec: float, econ_capex: float, econ_life: float, econ_dr: float,
             nu_mult: float, f_mult: float,
             base_env: dict) -> dict:
    ensure_dir(workdir)
    env = os.environ.copy()
    env.update(base_env)
    # Simulation horizon
    env['SIM_YEARS'] = str(years)
    env['TIME_STEP_S'] = str(timestep_s)
    # Weather
    if weather:
        env['LOAD_WEATHER_CSV'] = str(weather)
    # Economics
    env['ECON_ELEC_PRICE'] = str(econ_elec)
    env['ECON_CAPEX'] = str(econ_capex)
    env['ECON_LIFETIME_Y'] = str(econ_life)
    env['ECON_DR'] = str(econ_dr)
    # Enhancement bottom-only
    z_end = depth_m
    z_start = max(0.0, depth_m - max(0.0, L_enh_m))
    env['EHEP_ENABLE'] = '1'
    env['EHEP_NU_MULT'] = str(nu_mult)
    env['EHEP_F_MULT'] = str(f_mult)
    env['EHEP_Z_START_M'] = str(z_start)
    env['EHEP_Z_END_M'] = str(z_end)

    t0 = time.time()
    proc = subprocess.run([str(exe)], cwd=str(workdir), env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    elapsed_s = time.time() - t0
    out = proc.stdout
    # Paths
    results_csv = workdir / 'results.csv'
    summary_csv = workdir / 'summary.csv'

    # Parse summary
    summ = {}
    if summary_csv.exists():
        with summary_csv.open('r', encoding='utf-8', errors='ignore') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # only one row
                for k, v in row.items():
                    try:
                        summ[k] = float(v)
                    except Exception:
                        summ[k] = v
                break

    # Compute unmet from results
    unmet_kWh = 0.0
    if results_csv.exists():
        with results_csv.open('r', encoding='utf-8', errors='ignore') as f:
            reader = csv.DictReader(f)
            # derive dt_h from TIME_STEP_S
            dt_h = float(timestep_s) / 3600.0
            for row in reader:
                try:
                    unmet_kWh += float(row.get('Q_unmet_kW', '0') or 0) * dt_h
                except Exception:
                    pass

    # Objectives
    elec_kWh = (summ.get('comp_kWh', 0.0) or 0.0) + (summ.get('pump_kWh', 0.0) or 0.0)
    heat_kWh = (summ.get('heat_kWh', 0.0) or 0.0)
    scop_sys = heat_kWh / elec_kWh if elec_kWh > 1e-9 else 0.0
    res = {
        'L_enh_m': L_enh_m,
        'z_start_m': z_start,
        'z_end_m': z_end,
        'LCOH_annual': summ.get('LCOH_annual', 0.0) or 0.0,
        'elec_kWh': elec_kWh,
        'comp_kWh': summ.get('comp_kWh', 0.0) or 0.0,
        'pump_kWh': summ.get('pump_kWh', 0.0) or 0.0,
        'heat_kWh': heat_kWh,
        'SCOP_sys': scop_sys,
        'unmet_kWh': unmet_kWh,
        'elapsed_s': elapsed_s,
    }
    return res


def pareto_non_dominated(items, keys_min):
    # items: list of dicts; keys_min: objectives to minimize
    nd = []
    for i, a in enumerate(items):
        dominated = False
        for j, b in enumerate(items):
            if i == j:
                continue
            better_or_equal = True
            strictly_better = False
            for k in keys_min:
                if (b[k] > a[k]) and not math.isclose(b[k], a[k]):
                    better_or_equal = False
                    break
                if b[k] < a[k] and not math.isclose(b[k], a[k]):
                    strictly_better = True
            if better_or_equal and strictly_better:
                dominated = True
                break
        if not dominated:
            nd.append(a)
    return nd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exe', required=True, help='Path to enhance.exe')
    parser.add_argument('--weather', required=False, help='Path to yearly weather CSV (with T_out_C column)')
    parser.add_argument('--depth', type=float, help='Well depth (m). If omitted, read from DataConfig.h')
    parser.add_argument('--length-grid', default='0:800:9', help='Enhancement length grid: start:end:count or comma list (m)')
    parser.add_argument('--years', type=int, default=1)
    parser.add_argument('--dt', type=int, default=3600, help='Timestep in seconds (default 3600)')
    parser.add_argument('--econ-elec', type=float, default=0.0)
    parser.add_argument('--econ-capex', type=float, default=0.0)
    parser.add_argument('--econ-life', type=float, default=20.0)
    parser.add_argument('--econ-dr', type=float, default=0.08)
    parser.add_argument('--nu-mult', type=float, default=1.25)
    parser.add_argument('--f-mult', type=float, default=1.41)
    parser.add_argument('--runs', default='runs', help='Runs output directory')
    args = parser.parse_args()

    exe = Path(args.exe).resolve()
    proj_root = Path(__file__).resolve().parents[1]
    depth_m = args.depth if args.depth is not None else read_depth_from_header(proj_root)
    weather = Path(args.weather).resolve() if args.weather else None

    L_values = parse_grid(args.length_grid)

    runs_root = Path(args.runs).resolve()
    ensure_dir(runs_root)

    base_env = {
        'EHEP_ENABLE': '1',
    }

    results = []
    for idx, Lm in enumerate(L_values):
        wdir = runs_root / f'L{int(round(Lm))}m'
        res = run_case(exe, wdir, weather, args.years, args.dt,
                       depth_m, Lm,
                       args.econ_elec, args.econ_capex, args.econ_life, args.econ_dr,
                       args.nu_mult, args.f_mult,
                       base_env)
        results.append(res)
        print(f'Ran L_enh={Lm:.1f} m -> LCOH={res["LCOH_annual"]:.4g}, elec_kWh={res["elec_kWh"]:.1f}, unmet_kWh={res["unmet_kWh"]:.1f}')
        sys.stdout.flush()

    # Write results
    out_csv = runs_root / 'opt_results.csv'
    if results:
        with out_csv.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)

    # Pareto set (minimize these keys)
    keys = ['LCOH_annual', 'elec_kWh', 'unmet_kWh']
    pareto = pareto_non_dominated(results, keys)
    out_pareto = runs_root / 'opt_pareto.csv'
    if pareto:
        with out_pareto.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(pareto[0].keys()))
            writer.writeheader()
            writer.writerows(pareto)
    print(f'Wrote {out_csv} and {out_pareto} ({len(pareto)} Pareto points).')


if __name__ == '__main__':
    main()

