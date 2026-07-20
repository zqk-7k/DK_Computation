# GOAT 显式溶剂簇到 DFT opt/freq 的安全交接

本说明只覆盖“构象结果审核”和“生成 DFT 输入”。这些脚本不会自动批准生产构象，也不会启动 ORCA。

## 1. 为什么还有一个结果门

用户对 F3_019 等结构、pD 和显式水模型的批准，是对**计算假设**的批准。GOAT 完成以后，具体保留哪些构象是一个新的、依赖实际结果的决定。因此必须检查能量窗、连接关系、立体化学、质子化状态、反离子和显式溶剂簇是否完整，不能用先前的一次性批准代替。

门控规则如下：

- 每个 `figure_id × scenario_id` 独立处理，绝不混合 DMSO/D2O 或不同显式溶剂数目的构型。
- 同一 Figure 的清单内所有 GOAT starts 必须全部正常完成，且 state、当前输入、输出和 ensemble 的哈希一致。
- 草稿默认采用 12 kcal mol⁻¹ 的保守能量窗，以及固定原子顺序、全簇原子的 Kabsch RMSD 去重。
- 草稿生成器只能写 `DRAFT_NOT_APPROVED`，不能写生产批准。
- DFT 建立器只接受 `APPROVED_FOR_OPT_FREQ`，并要求顶层和每个保留构象的检查均为 `PASS`。
- 显式 DMSO/D2O 和反离子全部进入电子结构计算；它们只在 `nmr_atom_map_v2.csv` 中被标为不属于 analyte NMR，不会从坐标中删除。
- 氘按源 start 的 `isotope_map.csv` 逐原子写成 `M=2.01410177812`。

## 2. 生成只读审核草稿

PowerShell 中执行：

```powershell
$P = '${PACKAGE_ROOT}'
py "$P\scripts\draft_goat_selection_gate.py" --project "$P" --figure-id F3_001
```

当 52 个 Figure 的所有 GOAT starts 都完成后，可以省略 `--figure-id`，一次生成所有 Figure 的草稿。若任何 start 未完成、输入改变或哈希不一致，命令会在写文件前失败。

草稿位置：

```text
jobs\F3_xxx\goat_ensemble_selection\goat_ensemble_selection_DRAFT.json
jobs\F3_xxx\goat_ensemble_selection\draft_candidates\<scenario_id>\*.xyz
```

脚本不会创建 `goat_ensemble_selection_approved.json`。

## 3. 人工审核并批准

每个场景分别检查候选 XYZ：

1. analyte 的连接关系、立体化学和酰胺构型是否正确；
2. 盐型、质子化状态、反离子及显式溶剂数量是否与该 scenario 完全一致；
3. 坐标确实对应声明的 ensemble frame，且没有原子丢失或严重解离；
4. 能量窗和保留构象是否合理；
5. 删除被拒绝的候选，不要把拒绝项改成批准。

审核通过后，把草稿复制为：

```text
goat_ensemble_selection_approved.json
```

并且只在实际完成审核后填写：

```text
approval_status = APPROVED_FOR_OPT_FREQ
selection_completeness_check = PASS
protonation_and_cluster_integrity_review = PASS
reviewer = <审核人>
reviewed_at = <ISO 8601 时间>
```

每个保留构象也必须设置：

```text
approval_status = APPROVED_FOR_OPT_FREQ
connectivity_stereochemistry_check = PASS
cluster_composition_check = PASS
manual_extraction_check = PASS
```

不要改动候选 XYZ、GOAT ensemble、`atom_map.csv` 或 `isotope_map.csv`。任何字节变化都会使哈希验证失败。

## 4. 生成 DFT opt/freq 输入（仍不运行 ORCA）

```powershell
$P = '${PACKAGE_ROOT}'
py "$P\scripts\build_dft_jobs.py" --project "$P" --gate "$P\jobs\F3_001\goat_ensemble_selection\goat_ensemble_selection_approved.json"
```

如果不写 `--gate`，脚本会发现 `jobs\F3_*\...\goat_ensemble_selection_approved.json`。输出位于：

```text
02_dft_opt_freq_explicit\F3_xxx\<scenario_id>\conformer_xxx\
data\DFT_explicit_job_manifest.csv
```

每个 job 包含：

- `dft_opt_freq.inp`：`r2SCAN-3c CPCM(actual solvent) Opt TightOpt Freq TightSCF`；
- `approved_selected_cluster.xyz`：与审核候选字节一致；
- `source_atom_map.csv`、`source_isotope_map.csv`：来自该构象的原始 GOAT start；
- `nmr_atom_map_v2.csv`：analyte/environment 和有效同位素标记；
- `provenance.json`：gate、state、ensemble、坐标、map、manifest、scenario registry 和资源策略的哈希；
- `SCIENTIFIC_DEPENDENCY_SHA256.txt` 和明确的计算模型声明。

资源设置从 `config\resource_policy.json` 读取：本机为 6 MPI ranks、每 rank 4,000 MB，并保留最多 8 ranks 的硬上限。建立器不运行 ORCA；后续执行仍必须保持 `max-workers=1` 和一次只运行一个化合物/阶段。

## 5. 已有 ORCA 6.1.1 氘质量语法证据

只读测试使用：

```text
scratch\syntax_test_deuterium.inp
scratch\syntax_test_deuterium.out
```

已有输出正常结束，并报告总质量 20.03 AMU，证明本机 ORCA 6.1.1 接受 inline `M=2.01410177812`。本次交接测试只读取这两个文件，没有再次启动 ORCA。

## 6. 单元测试

```powershell
$P = '${PACKAGE_ROOT}'
py -m unittest discover -s "$P\tests" -v
```

测试完全使用临时的合成显式水簇，覆盖：DRAFT 不能进入 DFT、两个不同水合场景隔离、全簇原子保留、氘质量数目与 map 一致、DFT 输入 bundle 哈希与现有 runner 一致，以及读取既有 ORCA 氘语法证据。

