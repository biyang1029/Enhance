<#
  扫参脚本（PowerShell）— 放在项目根目录 enhance 下
  运行方法：
    1) 右键本文件 -> 使用 PowerShell 运行
    2) 或终端执行：powershell -NoProfile -ExecutionPolicy Bypass -File .\sweep.ps1
#>
Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

$ScriptDir = if ($PSScriptRoot) { $PSScriptRoot } else { Split-Path -Parent $MyInvocation.MyCommand.Path }

function Find-EnhanceExe([string]$baseDir){
  if ($env:ENHANCE_EXE -and (Test-Path $env:ENHANCE_EXE)) { return (Resolve-Path $env:ENHANCE_EXE).Path }
  $candidates = @(
    (Join-Path $baseDir 'x64\Release\enhance.exe'),
    (Join-Path $baseDir 'Release\enhance.exe'),
    (Join-Path $baseDir 'enhance.exe'),
    (Join-Path $baseDir 'x64\Debug\enhance.exe'),
    (Join-Path $baseDir 'Debug\enhance.exe')
  )
  foreach($p in $candidates){ if(Test-Path $p){ return (Resolve-Path $p).Path } }
  $found = Get-ChildItem -Path $baseDir -Filter 'enhance.exe' -Recurse -ErrorAction SilentlyContinue | Select-Object -First 1
  if($found){ return $found.FullName }
  throw '未找到 enhance.exe，请先编译 (Release|x64) 或设置 ENHANCE_EXE 环境变量。'
}

$exe = Find-EnhanceExe $ScriptDir
Write-Host "[Info] exe: $exe"

# 参数网格（可按需调整）
$hpList    = @(600)
$dtList    = @(0)           # 源侧每小时最大降温 (K/h)；0=关闭
$uaList    = @(22,30)
$mfList    = @(10)          # 源侧质量流量 (kg/s)
$effCarnotList    = @(0.65)
$evapApproachList = @(3.0)
$condApproachList = @(3.0)
$modRangeList     = @(1.0)

# 几何与材料
$dOuterList  = @(0.20)
$dInnerList  = @(0.10)
$boreList    = @(0.22)
$tPipeList   = @(0.008)
$soilKList   = @(2.5)
$soilRhoList = @(2600)
$soilCpList  = @(1600)
$groutKList  = @(4)
$pipeKInnerList = @(0.4)     # PE-RT II
$pipeKOuterList = @(4.0)     # J55 钢
$insulTopLenList = @(100)    # m 顶段等效保温长度
$insulKInnerList = @(0.1)
$enhBottomLenList= @(100)    # m 底部增强长度

$ts = Get-Date -Format 'yyyyMMdd_HHmmss'
$runRoot = Join-Path $ScriptDir ("runs_"+$ts)
New-Item -ItemType Directory -Force -Path $runRoot | Out-Null

# 基础环境
$env:SIM_YEARS   = '1'
$env:TIME_STEP_S = '600'
$env:LOG_HOURLY  = '1'

# 供暖季（10/15-4/15）
$env:HEAT_SEASON_ENABLE = '1'
$env:HEAT_START_MM = '10'
$env:HEAT_START_DD = '15'
$env:HEAT_END_MM   = '4'
$env:HEAT_END_DD   = '15'

# 天气与高温切断
$weather1 = Join-Path $ScriptDir 'weather_gansu.csv'
$weather2 = Join-Path (Join-Path $ScriptDir 'enhance') 'weather_gansu.csv'
if (Test-Path $weather1) { $env:LOAD_WEATHER_CSV = $weather1 }
elseif (Test-Path $weather2) { $env:LOAD_WEATHER_CSV = $weather2 }
else { $env:LOAD_WEATHER_CSV = 'weather_gansu.csv' }
$env:LOAD_CUTOFF_ENABLE = '1'
$env:LOAD_HEAT_CUTOFF_C = '26'
$env:LOAD_BASE_KW = '0'
$env:LOAD_INDOOR_T = '26'

# 热泵保护/调制
$env:HP_MIN_SRC_RETURN_C = '10'
$env:HP_MOD_RANGE_C      = '1'

# summary 头
$summaryFile = Join-Path $runRoot 'summary_runs.csv'
"hp_max_kW,UA_kWperK,flow_src_kgps,flow_load_kgps,set_load_out_C,tank_vol_m3,Q_load_kWh,Q_src_kWh,P_el_kWh,P_pump_kWh,dP_kPa_avg,avg_Q_load_kW,avg_Q_src_kW,COP_annual,COP_measured,COP_delta,HP_on_hours,HP_on_measured,HP_on_delta" | Set-Content -Path $summaryFile -Encoding UTF8

function Run-OneCombo([double]$hp,[double]$dt,[double]$ua,[double]$mf,[double]$effCar,[double]$evapK,[double]$condK,[double]$modC,[double]$pkIn,[double]$pkOut,[double]$insTop,[double]$insKin,[double]$enhB,[double]$dOuter,[double]$dInner,[double]$bore,[double]$tPipe,[double]$soilK,[double]$soilRho,[double]$soilCp,[double]$groutK){
  # 写环境变量
  $env:HP_MAX_Q_OUT_KW     = [string]$hp
  $env:HP_MAX_SRC_DT_PER_H = [string]$dt
  $env:LOAD_UA_KW_PER_K    = [string]$ua
  $env:MASS_FLOW_KGPS      = [string]$mf
  $env:HP_EFF_CARNOT       = [string]$effCar
  $env:HP_EVAP_APPROACH_K  = [string]$evapK
  $env:HP_COND_APPROACH_K  = [string]$condK
  $env:HP_MOD_RANGE_C      = [string]$modC
  $env:D_OUTER_M           = [string]$dOuter
  $env:D_INNER_M           = [string]$dInner
  $env:BORE_D_M            = [string]$bore
  $env:PIPE_THICK_M        = [string]$tPipe
  $env:SOIL_K              = [string]$soilK
  $env:SOIL_RHO            = [string]$soilRho
  $env:SOIL_CP             = [string]$soilCp
  $env:GROUT_K             = [string]$groutK
  $env:PIPE_K_INNER        = [string]$pkIn
  $env:PIPE_K_OUTER        = [string]$pkOut
  $env:INSUL_TOP_LEN_M     = [string]$insTop
  $env:INSUL_K_INNER       = [string]$insKin
  $env:EHEP_BOTTOM_LEN_M   = [string]$enhB

  # 参考值（实测）
  $HP_on_measured = 2890.0
  $COP_measured   = 4.3

  Write-Host "[Run] hp=$hp, dTcap=$dt, UA=$ua, mflow=$mf, kin=$pkIn, kout=$pkOut, ins=($insTop m,$insKin), enhB=$enhB"
  Push-Location $ScriptDir
    $env:TANK_VOL_M3 = '8'
    $env:TANK_SET_C  = '40'
    if (-not $env:FLOW_SRC_KGPS)  { $env:FLOW_SRC_KGPS  = [string]$mf }
    if (-not $env:FLOW_LOAD_KGPS) { $env:FLOW_LOAD_KGPS = '22.2' }
    & $exe
  Pop-Location
  if($LASTEXITCODE -ne 0){ Write-Warning "进程退出码 $LASTEXITCODE" }

  # 归档输出
  $tagDt = ([string]$dt) -replace '\.','p'
  $tag = "hp_${hp}_dt${tagDt}_ua${ua}_mf${mf}_ec${effCar}_ea${evapK}_ca${condK}_mr${modC}_kin${pkIn}_kout${pkOut}_ins${insTop}m_${insKin}_enh${enhB}m_do${dOuter}_di${dInner}_bd${bore}_tp${tPipe}"
  $resOut = Join-Path $runRoot ("results_"+$tag+".csv")
  $dbgOut = Join-Path $runRoot ("debug_"  +$tag+".csv")
  if(Test-Path (Join-Path $ScriptDir 'results.csv')){ Copy-Item (Join-Path $ScriptDir 'results.csv') $resOut -Force }
  if(Test-Path (Join-Path $ScriptDir 'debug.csv'))  { Copy-Item (Join-Path $ScriptDir 'debug.csv')  $dbgOut -Force }

  # 汇总（每小时一行，直接累加即为 kWh）
  if(Test-Path $resOut){
    $rows = Import-Csv -Path $resOut
    if($rows.Count -gt 0){
      $sumQspace=0.0; $sumQdhw=0.0; $sumPel=0.0; $sumPpump=0.0; $sumQsrc=0.0
      $sumDP=0.0; $sumQoutInst=0.0; $sumQsrcInst=0.0; $onRows=0
      $sumFlowSrc=0.0; $sumFlowLoad=0.0
      foreach($r in $rows){
        $qs = [double]$r.Q_space_served_kW
        $qd = [double]$r.Q_dhw_served_kW
        $pe = [double]$r.P_el_kW
        $pp = [double]$r.P_pump_kW
        $qg = 0.0; if($r.PSObject.Properties.Name -contains 'Q_geo_kW'){ $qg = [double]$r.Q_geo_kW }
        $sumQspace += $qs
        $sumQdhw   += $qd
        $sumPel    += $pe
        $sumPpump  += $pp
        $sumQsrc   += $qg
        $on = 0; try { $on = [int][double]$r.HP_on } catch {}
        if($on -gt 0){
          $onRows += 1
          if($r.PSObject.Properties.Name -contains 'dP_kPa'){ $sumDP += [double]$r.dP_kPa }
          if($r.PSObject.Properties.Name -contains 'Q_out_kW'){ $sumQoutInst += [double]$r.Q_out_kW }
          if($r.PSObject.Properties.Name -contains 'Q_geo_kW'){ $sumQsrcInst += [double]$r.Q_geo_kW }
          if($r.PSObject.Properties.Name -contains 'flow_src_kgps'){ $sumFlowSrc += [double]$r.flow_src_kgps }
          if($r.PSObject.Properties.Name -contains 'flow_load_kgps'){ $sumFlowLoad += [double]$r.flow_load_kgps }
        }
      }
      $Q_load_kWh = $sumQspace + $sumQdhw
      $Q_src_kWh  = $sumQsrc
      $P_el_kWh   = $sumPel
      $P_pump_kWh = $sumPpump
      $COP_annual = if(($P_el_kWh + $P_pump_kWh) -gt 1e-9){ [double]($Q_load_kWh / ($P_el_kWh + $P_pump_kWh)) } else { 0.0 }
      $HP_on_hours = [double]$onRows
      $dP_avg = if($onRows -gt 0){ [double]($sumDP / $onRows) } else { 0.0 }
      $avg_Q_load_kW = if($onRows -gt 0){ [double]($sumQoutInst / $onRows) } else { 0.0 }
      $avg_Q_src_kW  = if($onRows -gt 0){ [double]($sumQsrcInst / $onRows) } else { 0.0 }
      $flow_src_kgps = if($onRows -gt 0){ [double]($sumFlowSrc / $onRows) } else { [double]$mf }
      $flow_load_kgps= if($onRows -gt 0){ [double]($sumFlowLoad / $onRows) } else { 22.2 }
      $COP_delta = [double]($COP_annual - $COP_measured)
      $HP_on_delta = [double]($HP_on_hours - $HP_on_measured)
      Add-Content -Path $summaryFile -Value ("$hp,$ua,$flow_src_kgps,$flow_load_kgps,40,8,$Q_load_kWh,$Q_src_kWh,$P_el_kWh,$P_pump_kWh,$dP_avg,$avg_Q_load_kW,$avg_Q_src_kW,$COP_annual,$COP_measured,$COP_delta,$HP_on_hours,$HP_on_measured,$HP_on_delta")
    }
  }
}

foreach($hp in $hpList){
  foreach($dt in $dtList){
    foreach($ua in $uaList){
      foreach($mf in $mfList){
        foreach($kin in $pipeKInnerList){
          foreach($kout in $pipeKOuterList){
            foreach($itop in $insulTopLenList){
              foreach($ikin in $insulKInnerList){
                foreach($enhB in $enhBottomLenList){
                  foreach($do in $dOuterList){
                    foreach($di in $dInnerList){
                      foreach($bd in $boreList){
                        foreach($tp in $tPipeList){
                          foreach($sk in $soilKList){
                            foreach($sr in $soilRhoList){
                              foreach($scp in $soilCpList){
                                foreach($gk in $groutKList){
                                  foreach($ec in $effCarnotList){
                                    foreach($ea in $evapApproachList){
                                      foreach($ca in $condApproachList){
                                        foreach($mr in $modRangeList){
                                          Run-OneCombo -hp $hp -dt $dt -ua $ua -mf $mf -effCar $ec -evapK $ea -condK $ca -modC $mr -pkIn $kin -pkOut $kout -insTop $itop -insKin $ikin -enhB $enhB -dOuter $do -dInner $di -bore $bd -tPipe $tp -soilK $sk -soilRho $sr -soilCp $scp -groutK $gk
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

Write-Host "[Done] 输出目录: $runRoot"

