# opt/freq、QRRHO 权重与 NMR 属性层的安全流程

本层仍然采用“默认只检查、不运行”的原则。新增脚本不会自动批准任何生产结果。本次实现和测试没有启动生产 ORCA。

## 1. opt/freq 串行 runner

先做无计算的完整哈希检查：

```powershell
$P = '${PACKAGE_ROOT}'
py "$P\scripts\run_opt_freq_serial.py" --project "$P"
```

默认输出应为 `DRY RUN`，不会启动 ORCA。它会重新核对：批准的 GOAT 构象门、GOAT state/ensemble、选中坐标、origin atom/isotope maps、DFT input、科学依赖哈希和当前 manifest。

只有在构象门已经人工批准且准备开始单个任务时，才使用：

```powershell
py "$P\scripts\run_opt_freq_serial.py" --project "$P" --only F3_xxx__scenario__c001 --execute
```

runner 具有以下硬约束：

- 全局独占锁；
- 一次只运行一个 job；
- 运行前检查外部 ORCA/MPI/xTB 进程；
- 按 `resource_policy.json` 采样 CPU，并检查 RAM 和磁盘；
- opt/freq 使用 6 MPI ranks、每 rank 4,000 MB，绝对上限 8 ranks；
- 运行时持续 CPU 高于 88% 达 90 s 时，完成当前 job 后停止，不强制杀死健康 ORCA；
- 每个 attempt 的 input、output、GBW、Hessian、property 和 optimized XYZ 均记录 SHA-256。

生产执行命令在本次工作中**没有调用**。

## 2. 生成 minimum/QRRHO/Boltzmann 草稿

所有 opt/freq state 完成后：

```powershell
py "$P\scripts\draft_minimum_boltzmann_gate.py" --project "$P" --figure-id F3_xxx
```

脚本读取 ORCA 6.1 输出中明确声明的 Grimme QRRHO 热化学结果，要求 100.0 cm⁻¹ reference frequency，并提取：

- 最终电子能；
- ORCA 报告的 QRRHO Gibbs free energy；
- 热化学温度；
- 全部频率和虚频数目；
- optimized XYZ、Hessian、atom map 和 state 哈希。

权重按同一 `figure_id × scenario_id × nucleus × temperature` 归一化：

```text
w_i = exp[-(G_i - G_min)/(RT)] / Σ_j exp[-(G_j - G_min)/(RT)]
```

不同溶剂、显式溶剂数、质子化态和反离子构型不会混合。当前 DFT thermochemistry 锁定在 298.15 K；如果真实采集温度不同，必须提供匹配温度的热化学结果，不能只改 JSON 中的温度。

输出是：

```text
jobs\F3_xxx\minimum_boltzmann_audit\minimum_boltzmann_DRAFT.json
```

它永远是 `DRAFT_NOT_APPROVED`。

## 3. 人工审核 minimum gate

逐个 opt/freq state 检查：

1. 零虚频；
2. 未发生错误质子转移；
3. analyte 连接、立体化学和盐型未改变；
4. 显式溶剂和反离子没有不合理逃逸；
5. 每个 scenario 至少保留一个构象；
6. 删除构象后，重新按该 scenario/nucleus/temperature 计算并归一化权重。

批准文件名为 `minimum_boltzmann_approved.json`。顶层使用 `APPROVED_FOR_NMR`，`same_composition_only_check=PASS`；保留构象也使用 `APPROVED_FOR_NMR` 且两个人工检查为 `PASS`。拒绝项使用 `REJECTED_AFTER_OPT_FREQ` 并写明理由。每个 weight set 也必须为 `APPROVED_FOR_NMR`。

## 4. 当前仍然阻止属性计算的三个显式 gate

以下文件目前有意保持未批准：

- `config\reference_job_registry.json`：reference input/output、方法、溶剂域和实验 reference convention 尚未锁定；
- `jobs\j_method_benchmark_gate.json`：PBE/pcJ-3 尚未在代表性 H–H、H–F/F–C 和 H–P/P–C 集合上与 PBE0/pcJ-2 比较；
- `jobs\F3_006_iodine_x2c_gate.json`：碘的 basis、picture-change、X2C sensitivity 和 J basis 兼容性尚未通过。

`build_property_jobs.py` 遇到任何一个相关的 DRAFT gate 都会在写生产输入前失败。不能用前面的结构/场景总体批准代替这些结果依赖和方法依赖的批准。

## 5. reference registry 的严格匹配

每个 ¹H/¹³C/¹⁹F/³¹P reference 必须同时匹配：

- nucleus；
- `DMSO-d6` 或 `D2O` solvent domain；
- 完整 shielding method fingerprint；
- reference input/output 文件及 SHA-256；
- 实际实验 reference convention。

只要 method fingerprint 或任一哈希不一致，属性 builder 就会拒绝。不同溶剂域的 reference 不能共用。

## 6. 构建 analyte-only 属性输入

上述 gate 均真实批准以后：

```powershell
py "$P\scripts\build_property_jobs.py" --project "$P" --gate "$P\jobs\F3_xxx\minimum_boltzmann_audit\minimum_boltzmann_approved.json"
```

脚本只生成输入，不启动 ORCA。输出包括：

```text
03_nmr_properties_explicit\F3_xxx\<scenario>\<conformer>\shielding_*\
03_nmr_properties_explicit\F3_xxx\<scenario>\<conformer>\j_coupling\
data\NMR_shielding_job_manifest.csv
data\NMR_j_coupling_job_manifest.csv
data\NMR_conformer_weights.csv
```

电子结构模型仍包含整个 analyte + 显式溶剂 + 反离子簇。`nmr_atom_map_v2.csv` 中不属于 analyte 的核不会出现在 shielding 或 J 请求中。

ORCA 6.1 的索引规则并不统一：

- `Nuclei = ...` 从 1 开始；
- `SpinSpinAtomPairs {...}` 从 0 开始。

因此每个 J job 都附有 `j_pair_allowlist.csv`，同时保存 1-based 和 0-based 编号。实现依据是 [ORCA 6.1 官方 NMR 手册](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/nmr.html)。

J allowlist 只保留 analyte 内、6.0 Å 以内的 H–H、H–F、F–F、C–F、H–P、P–P 和 C–P 类型；环境 DMSO/D2O 和反离子核不会进入 perturbation list。

## 7. 测试

```powershell
py -m unittest discover -s "$P\tests" -v
```

合成夹具覆盖：单实例锁、依赖篡改拒绝、场景不混合、QRRHO 权重逐组归一、environment nuclei 排除、ORCA 两套索引、reference method mismatch、DRAFT reference registry 和 DRAFT J benchmark fail-closed。测试不运行 ORCA。

