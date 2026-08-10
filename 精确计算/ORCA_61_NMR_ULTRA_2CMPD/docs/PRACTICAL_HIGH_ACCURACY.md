# PRACTICAL_HIGH_ACCURACY：默认含痕量水的生产模式

> 状态：`INPUTS_PREPARED / HPC_NOT_RUN / DRAFT_NOT_VALIDATED`

本模式用于日常 `1H`/`13C` NMR 预测。`ULTRA_DIAGNOSTIC` 保留为异常结果的升级诊断模式。两种模式都只生成计算结果，不能替代实测表征。

## 1. 默认计算范围

苯甲酰胺使用 3 个、`(S)`-脯氨酰胺使用 4 个相互隔离的固定计量分支，每个分支 10 个确定性 GOAT 起点。第二个分子多出的 `2 DMSO` 无水分支用于和 `2 DMSO + 1 H2O` 做同组成水效应比较：

| 化合物 | 主连续介质 | 无水第一溶剂壳 | 含水分支的同组成无水对照 | 默认必算的含水分支 |
|---|---|---|---|---|
| 苯甲酰胺 `F3_005` | `DMSO_00_W00_continuum` | `DMSO_02_W00_primary` | 同左 | `DMSO_02_W01_trace_water_sensitivity` |
| `(S)`-脯氨酰胺 `F3_015` | `DMSO_00_W00_continuum` | `DMSO_03_W00_extended` | `DMSO_02_W00_primary` | `DMSO_02_W01_trace_water_sensitivity` |

默认 practical manifest 共 70 个 GOAT 作业，其中 20 个属于 `2 DMSO + 1 H2O`。含水分支影响构象、能量和分析物屏蔽，但显式水原子保持 `environment_atom=true`，不会作为分析物核输出。

含水和无水分支必须分别完成构象选择、opt/freq、QRRHO 权重和 NMR 性质计算。水效应只比较 `DMSO_02_W01` 与同一化合物的 `DMSO_02_W00`，避免把水效应与 DMSO 数量效应混淆。禁止把不同显式溶剂计量直接做普通 Boltzmann 混合。含水分支的角色是条件敏感性区间，不是实测水活度或主结果的自动替代。

两个匹配分支完成后，生成逐原子的水敏感性差值表：

```bash
python3 scripts/compare_trace_water_sensitivity.py \
  --dry-bundle NMR_predictions_CALCULATED_NOT_EXPERIMENTAL/F3_005/DMSO_02_W00_primary/F3_005__DMSO_02_W00_primary_prediction_bundle.json \
  --wet-bundle NMR_predictions_CALCULATED_NOT_EXPERIMENTAL/F3_005/DMSO_02_W01_trace_water_sensitivity/F3_005__DMSO_02_W01_trace_water_sensitivity_prediction_bundle.json \
  --output-directory deliverables/F3_005/trace_water_sensitivity
```

该表只报告 `wet − dry` 的 Δδ 和环境敏感性标记。比较前会强制核对分子身份、opt/freq–QRRHO 协议指纹、shielding 方法、TMS/参考注册表、核种温度、权重政策和溶液等价化配置。非交换 `1H`、`13C` 和交换 `1H` 的默认触发阈值分别是 0.10、1.0 和 0.50 ppm；交换氢次序改变也会触发升级。脚本禁止使用不匹配的干/湿组成，也不会根据实验吻合度替换主结果。

生成/刷新两个可运行的 profile manifest：

```bash
python3 scripts/build_workflow_profile.py
python3 scripts/build_workflow_profile.py --profile ULTRA_DIAGNOSTIC
```

默认先以 practical 做预检并渲染小规模 pilot：

```bash
PROFILE_MANIFEST=data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/GOAT_explicit_job_manifest.csv

python3 hpc/bin/preflight.py \
  --package-root "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --stage goat \
  --scheduler slurm

hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 1-2 --max-concurrent 2
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 21-22 --max-concurrent 2
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 51-52 --max-concurrent 2
```

三个命令分别覆盖普通无水 pilot，以及两个化合物各自的含水 pilot；默认都只渲染调度脚本，人工检查后才追加 `--submit`。含水分支需要检查水分子身份、团簇完整性、酰胺/水质子转移、SCF 尾部和虚频。若优化后发生质子化状态改变，必须建立新物种，不能继续当作普通构象。

pilot 通过后，初次全量建议 `--max-concurrent 16`。GOAT 每个任务申请 8 CPU 和 24 GiB；只有当实测内存、scratch、许可证和队列均允许时才逐步升到 32 或 60。60 并发需要 480 CPU，并应至少有约 1.44 TiB 可调度内存。并发上限按 `min(floor(可用CPU/8), floor(可用内存/24 GiB), 站点/许可证上限, scratch实测上限, 70)` 计算；增加并发只缩短批次墙钟时间，不提高单任务精度。

## 2. 自适应性质层

- `PBE0/pcSseg-3`：全部批准构象；
- `revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore`：只覆盖 `1H` 与 `13C` 各自权重累计 95% 所需构象的并集；
- J：只在 J 方法门控通过后运行，并使用 practical 人口范围；
- RI-MP2、DLPNO-MP2、AutoCI-CCSD(T) 和 VPT2：仍为人工触发的诊断分支。

构建性质作业时必须显式绑定 profile：

```bash
python3 scripts/build_property_jobs.py \
  --workflow-profiles config/workflow_profiles.json \
  --workflow-profile PRACTICAL_HIGH_ACCURACY
```

## 3. 三层谱图交付

每个 `figure_id/scenario_id/nucleus/shift_variant` 必须分别保留：

1. `CLEAN_CALCULATED`：仅计算峰位、计算 J 和显示用干净线型；
2. `CONDITIONED_DETERMINISTIC`：经验线宽、DMSO/水峰和确定性基线；
3. `CONDITIONED_SEEDED_REALIZATION`：在第 2 层上加入可复现的随机噪声。

制谱命令会先输出规范 `*_appearance_lines.csv`，再调用独立外观层：

```bash
python3 scripts/generate_nmr_spectra.py \
  --render-appearance \
  --appearance-config config/appearance_profiles.json
```

外观层不得移动任何分析物 ppm，也不得再次分裂已经由自旋 Hamiltonian 生成的 resolved transitions。输出路径包含 `scenario_id`，所以含水和无水结果不会互相覆盖。

## 4. 苯甲酰胺实测条件化外观

`F3_005` 默认使用两个单样本条件化 profile：

- `benzamide_1_1_dmso_d6_1H_400MHz_calibrated_single_sample`；
- `benzamide_1_1_dmso_d6_13C_100MHz_calibrated_single_sample`。

它们只从 `real_spectrum_calibration/benzamide_1-1_calibration.json` 读取仪器/处理条件、噪声尺度、基线尺度、线宽以及 DMSO/水峰外观。配置锁定该文件 SHA-256，并明确禁止逐峰拟合计算位移。

本地 Bruker raw FID 后来已单独核验并在 `validation/benzamide_raw_fid_provenance.json` 记录哈希，但未随公开包分发。当前外观参数仍来自处理后的 PDF/XY 与嵌入采集元数据，没有从 receiver-domain FID 重新拟合；元数据会把“外观参数校准”和“分析物位移拟合”分开记录。

已条件化的主要值包括：

- `1H`：400.132 MHz、312.609 K、非交换峰 FWHM 3.33 Hz、交换峰总 FWHM 18.16 Hz、DMSO 2.5004 ppm、水 3.3172 ppm；
- `13C`：100.613 MHz、304.065 K、分析物默认 FWHM 2.5 Hz、DMSO-d6 七重峰中心 39.5199 ppm、间隔 21.1 Hz、分峰 FWHM 6.43 Hz；
- 噪声与基线以分析物最大峰为参考，不随水峰或溶剂峰高度被动放大；
- pseudo-Voigt 数值实现保证配置 FWHM 就是最终曲线 FWHM。

这些参数只适用于苯甲酰胺样品 `1-1` 的条件化外观，不能静默迁移到其他药物。`F3_015` 仍使用通用 provisional profile，直到取得它自己的采集元数据或可审计数值谱。

## 5. 科学边界

显式水微态与外观水峰是两个不同概念：

- 显式 `H2O` 是静态结构/屏蔽敏感性模型；
- 3.3172 ppm 水峰是经验采集外观；
- 外观水峰高度不能决定湿微态人口；
- 静态 ORCA 不能预测噪声、基线、残余溶剂强度、真实峰高、弛豫线宽或交换展宽。

所有外观图必须显示 `SYNTHETIC APPEARANCE - NOT EXPERIMENTAL`，并与纯计算谱同时交付。
