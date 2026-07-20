# ORCA 52 种化合物高质量串行流程：如何调用

状态：`LOCAL RUNNER VERIFIED / INPUT JOBS NOT YET APPROVED / NO NEW ORCA CALCULATION STARTED`

这不是一段抽象提示词，而是一套可从本机 PowerShell 直接调用的安全串行框架。默认只预览；只有显式加入 `-Execute`，并且结构、物种、pD、采集条件、reference、3D 坐标、atom/isotope map、方法、哈希和阶段顺序全部通过审核，才允许启动 ORCA。

## 1. 已锁定的溶剂与主要物种

- **F3_019（L-asparagine）**：D2O 主条件，L-两性离子为主模型。实际 pD 优先；若确实没有记录，可用 pD 7.4/DSS 作为标准化假设，但必须明确标为 `ASSUMED`。水合数测试 6、8、10。
- **F3_020（L-glutamine）**：D2O 主条件，近中性 pD 下的 L-两性离子为主模型。图中的中性式只做独立物种敏感性计算，不与两性离子做 Boltzmann 混合。水合数测试 6、8、10。
- **F3_021（L-argininamide dihydrochloride）**：D2O 主条件。酸性条件主模型为二阳离子 + 2 个溶剂分隔的 Cl−；接触离子对及较高 pD 下的单阳离子各自单独计算。无实际 pD 时，可用 pD 5.5 作为二盐酸盐假设并清楚标注。水合数测试 8、10、12。
- **F3_035**：D2O。近中性 pD 主模型为 sulfamate anion；普通 ^1H 只有 HOD/背景层，不把交换性 NH 当作稳定产物峰；无 ^13C。水合数测试 4、6、8。
- **F3_043**：D2O。近 pD 7 以中性 phosphoramide 为主，^31P{^1H} 是主要鉴定核；^1H 为背景/交换控制，无 ^13C。必须记录溶液放置时间并检查水解。
- **F3_044**：D2O。近 pD 7 以 monoanionic diamidophosphate 为主，记录真实反离子；^31P 为主要鉴定核，^1H 为背景/交换控制，无 ^13C。必须做时间零点/延迟 ^31P 水解检查。

F3_019/020/021 若要补 DMSO-d6，只能建立另一个 condition、manifest、构象系综、reference 和输出命名空间（如 `F3_019__DMSO`）。D2O 与 DMSO-d6 不是同一系综，绝不混合平均。

详细规则见：

- `config/condition_overrides.json`
- `..\ORCA_52_NMR_protocol_for_review\Figure3_52_proposed_NMR_conditions.csv`
- `..\ORCA_52_NMR_protocol_for_review\ORCA_52_NMR_master_instruction.md`

## 2. 真实谱图校准结论

已只读盘点 `${LOCAL_SOURCE_NOT_INCLUDED}\处理过的谱图`：27 个 Mnova 容器、27 个 PDF、26 个 TIFF 和 1 个处理后 XY CSV。可恢复的 29 条采集/处理记录中，24 条为 DMSO、5 条为 CDCl3；16 条 ^1H 均为 NS=8、D1=1 s、LB=0.3 Hz、SI=65536，13 条 ^13C 均为 NS=512、D1=2 s、LB=1.0 Hz、SI=32768。

苯甲酰胺 `1-1-H/C` 的本地校准包括：

- ^1H：400.132470802 MHz，312.609 K；残余 DMSO 2.5004 ppm，水 3.3172 ppm；非交换峰约 3.3 Hz，DMSO 未分辨包络约 6.7 Hz，交换峰可约 18.2 Hz。
- ^13C：100.622829803 MHz，304.0651 K；产物峰 FWHM 约 1.77–4.45 Hz；DMSO-d6 七重线间隔约 21.07 Hz，单线宽约 6.31–6.54 Hz。
- ^1H 与 ^13C 使用不同噪声、线宽、基线、动态范围和溶剂峰模型。每张合成外观谱使用独立、可复现的噪声随机种子。

当前文件夹中没有可独立核验的 Bruker `fid/ser/acqus/procs` 原始目录。`.mnova`、PDF、TIFF 和处理后 CSV 不能自动称为 vendor FID；少数 Mnova 元数据保存了以 `/fid` 结尾的历史路径，也不等于原始时域数据仍在。

完整盘点见：

- `..\ORCA_52_NMR_protocol_for_review\real_spectrum_calibration\processed_spectra_inventory.csv`
- `..\ORCA_52_NMR_protocol_for_review\real_spectrum_calibration\benzamide_calibration_notes.md`
- `real_spectrum_calibration\benzamide_1-1_calibration.json`

若 Mnova 已保存并关闭，可运行下面命令逐个检查容器能否通过 Mnova API 导出内部时域数组和处理后 XY；脚本不会修改源文件：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\scripts\Export-MnovaCalibration.ps1"
```

只有 `fidCount()>0` 且实际导出复数时域数组，才能写成“Mnova 容器中可导出的内部 FID”；仍需与原始 Bruker 实验目录核对，才能称为原始 vendor FID。

## 3. 本机资源策略

本机已审计为 Lenovo 82TD、i9-12900HX（16 物理核/24 逻辑线程）、63.73 GiB RAM。ORCA 6.1.1 位于 `${ORCA_EXE}`；xTB 6.7.1pre 位于 `${XTB_DIR}\xtb.exe`；MS-MPI 为 10.1.12498.18。三者路径、版本与 SHA-256 已写入 `config/resource_policy.json`，每次启动前复核。

固定策略：

- 同时只计算 1 个化合物、1 个阶段、1 个构象；
- 每个固定团簇计量均从至少 10 个不同、哈希不相同的起始 3D 构象/溶剂排布逐个做 GOAT；输入固定 `RandomSeed false` 以便确定性复现，但不把这些搜索误称为编号的 GOAT 随机种子；
- GOAT 为 4 MPI ranks、`NWorkers 4`、`%maxcore 1200`；每个固定 scenario/显式溶剂计量的 10 个起点也全部逐个运行；
- 普通 opt/freq、shielding、J 分别最多 6 MPI ranks，OMP/MKL/OpenBLAS 均为 1 线程；I/X2C 或大水合团簇用 4 ranks；
- 启动前对 CPU 采样 30 s；平均负载高于 45% 不启动；CPU 高于 88% 持续 90 s 时，安全完成当前 ORCA 阶段后停止，不强杀健康任务；
- 这是“低 ranks + 单实例 + 软停止”策略，不是操作系统级硬性 CPU 限速。正常 4–6/24 ranks 通常远低于 90%，但仍以实时状态为准；
- C/D/E 剩余空间低于 20/20/100 GiB 时不启动；每次 attempt 使用独立 scratch，失败数据保留诊断；
- runner、输入、XYZ、atom map、isotope map、结构验证报告、condition/job lock 和 reference 都以 SHA-256 锁定；运行快照保存在 attempt 内。

## 4. 最直接的本地调用方式

PowerShell 执行策略会阻止直接双击 `.ps1`，请复制完整命令。

### 4.1 查看状态（不计算）

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_serial.ps1" -Status
```

### 4.2 预览一个化合物（不计算）

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_serial.ps1" -Compound F3_019
```

可把 `F3_019` 换成任意编号。F3_005 是已经完成的苯甲酰胺 demo，在队列中禁用并跳过。

### 4.3 审核完成后，只运行一个化合物

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_serial.ps1" -Compound F3_019 -Execute -Resume
```

runner 现在支持分阶段 manifest。第一阶段只需准备 `condition_lock → 至少10个GOAT起点 → goat_ensemble_selection`；此时不要求依赖GOAT结果的 opt/freq、shielding 或 J 输入预先存在。也可以在完整 manifest 中把下游 ORCA 行标为 `PENDING_UPSTREAM`。当运行真正到达该行时，runner 才检查它的输入、坐标、job lock、上游人工 gate 及全部 SHA-256；未准备完成会安全停止且不会启动ORCA。

第一阶段可直接从 `jobs\stage_manifest_phase1_goat_template.csv` 复制；完整累计模板为 `jobs\stage_manifest_template.csv`。分阶段manifest必须是“累计式”的：第二阶段是在第一阶段全部行之后追加opt/freq和`minimum_boltzmann_audit`，不能删除已经产生状态的上游行。

用户对整体计算计划的批准不能代替结果依赖型审核。GOAT完成后仍必须根据真实输出填写并批准 `goat_ensemble_selection_approved.json`；opt/freq完成后仍必须根据频率、团簇完整性、能量和温度填写并批准 `minimum_boltzmann_approved.json`。模板分别为：

- `jobs\goat_ensemble_selection_gate_template.json`
- `jobs\minimum_boltzmann_gate_template.json`

严格依赖链为：

`GOAT（每个scenario至少10个独立起点） → goat_ensemble_selection人工gate → opt/freq → minimum_boltzmann_audit人工gate → shielding/J → spin simulation → experimental appearance → QA`

GOAT阶段允许reference job仍为 `DRAFT_NOT_APPROVED`；但到达shielding阶段时，method-matched reference必须已经是`APPROVED`，且其输入、输出和方法哈希全部通过复核。

### 4.4 请求“当前阶段完成后停止”

在**第二个 PowerShell 窗口**运行：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_serial.ps1" -StopAfterCurrentStage
```

不要用关闭 runner 窗口或 Ctrl+C 代替正常停止；ORCA 子进程可能仍在后台完成，而 runner.lock 会被释放。即使如此，下一次启动的 external-process gate 仍会检测并拒绝重叠任务。

### 4.5 只验证分阶段门控合同（不会启动 ORCA）

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\tests\Test-StagedManualGateContract.ps1"
```

该测试解析 `run_serial.ps1` 的 PowerShell AST，检查第一阶段在 `goat_ensemble_selection` 结束、完整模板中的两个结果依赖型人工 gate 顺序、下游 `PENDING_UPSTREAM` 状态、opt/freq 与 shielding/J 的几何哈希来源，以及拒绝分支位于任何 `Invoke-OrcaStage` 调用之前。加 `-RuntimeSmoke` 只会附加运行 `-Status` 与 dry-run，仍不会启动 ORCA。

## 5. 每个化合物的逐项顺序

1. 从 Figure 3 registry 锁定母体结构、手性和 parent identity。
2. 按实验溶剂/pD 建立 scenario-specific 溶液物种；盐、质子化态、反离子和显式溶剂计量分别建模。
3. 分核记录实际或明确假设的温度、频率、reference 和采集参数；苯甲酰胺已经证明 ^1H/^13C 温度不同，不能用一个全局温度覆盖。
4. 锁定计划使用的reference ID；reference计算可以与GOAT/opt-freq分阶段完成，但在任何shielding任务启动前必须批准并锁定输入/输出 SHA-256。
5. 每个固定团簇计量生成至少 10 个独立起始 3D 构象/溶剂排布，逐个通过连接、手性、质子化、元素计数、charge/multiplicity、atom map 和 isotope map 验证。
6. 逐个GOAT完成后先暂停；人工从哈希锁定的ensemble中提取并批准构象，opt/freq输入必须逐字节使用这些批准的XYZ。
7. 逐个做r2SCAN-3c opt/freq，然后审核虚频、质子转移、团簇逃逸和离子对；仅零虚频且结构完整的极小值可进入按各核采集温度计算的QRRHO/Boltzmann权重。shielding/J输入必须逐字节使用该gate批准的优化XYZ。
8. 逐构象做 shielding；逐构象做需要的 J。F3_006 使用专用 I mixed-basis 与 X2C 敏感性分支。
9. 单线程完成 spin simulation；保存 stick、干净展宽谱和完整 J/shift 表。
10. 单线程生成实验外观层：加入与真实分布匹配但每个样品独立的噪声、基线、相位/数字伪影、DMSO/HOD/水背景和核种特异线宽。
11. 保存纯计算、背景分量和合成外观三套结果以及 TXT/CSV/JSON/SHA-256；通过最终 QA 后，才进入下一个化合物。

构象审核、spin simulation、实验外观和 QA 当前是强制人工 gate；专用 Python 后处理器尚未被冒充为“已经全自动完成”。

## 6. 科学边界

ORCA 计算的是屏蔽常数、J 耦合和相关量；HRMS 的 m/z/同位素包络与碎裂模型不是由 NMR 计算自动产生。最终输出必须永久标记：

`CALCULATED / SYNTHETIC APPEARANCE — NOT EXPERIMENTALLY ACQUIRED`

这些结果可用于谱峰归属、实验设计和方法评估，但不能代替向 Nature 投稿所需的真实原始表征数据。
