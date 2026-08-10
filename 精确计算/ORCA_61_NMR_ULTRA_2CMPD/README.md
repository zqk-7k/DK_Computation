# ORCA 6.1 高精度 NMR 双化合物工作流

> **当前状态：`INPUTS_PREPARED / HPC_NOT_RUN / DRAFT_NOT_VALIDATED`**
> 本目录已经准备了可审计的 ORCA 6.1.1 输入、门控规则和 Linux HPC 调度层，但尚未在真实超算上完成 ORCA 计算，也尚未通过独立实验验证。因此，当前包不能称为“已经验证的参数包”，目录中也没有可冒充实测数据的最终预测谱。

本包用于计算两个中性分子的 `¹H` 和 `¹³C` NMR：

| ID | 化合物 | 计算物种 | 分子式 | 主要条件 |
|---|---|---|---|---|
| `F3_005` | 苯甲酰胺 / benzamide | neutral | C₇H₇NO | DMSO-d₆ |
| `F3_015` | `(S)`-脯氨酰胺 / neutral `(S)`-prolinamide | neutral free base, `NC(=O)[C@@H]1CCCN1` | C₅H₁₀N₂O | DMSO-d₆ |

`F3_015` 按图示中性游离碱处理，不代表盐酸盐或其他盐。若实验样品实际为盐、不同质子化态或不同对映体，必须先更新结构与条件锁，不能直接使用当前结果。

## 设计目标与诚实边界

本包把“电子结构预测”和“仪器外观合成”分开：

1. ORCA 层预测构象、自由能、屏蔽常数和 J 耦合；
2. 参考层使用同方法 TMS 将屏蔽常数转换为化学位移；
3. 构象层只在**相同显式溶剂计量**内进行 Boltzmann 平均；
4. 外观层可加入 DMSO、水、基线、噪声、线宽和强度模型，但必须标记为 `SYNTHETIC / NOT EXPERIMENTAL`。

目标是尽可能接近相同条件下的实验峰位、顺序和多重性，而不是复制实验曲线。噪声、基线、溶剂峰、峰高、交换峰宽和数字处理历史不能由静态第一性原理唯一确定。任何外观层参数都不能用于反向移动分析物峰位，也不能把合成谱称为实测谱。

## 理论层级

生产候选主线为：

- 构象搜索：`GOAT-ENTROPY/XTB/ALPB(DMSO)`，每个固定计量场景 20 个确定性起点；
- 最终优化与频率：`PBE0-D4/def2-TZVPP/RIJCOSX/CPCM(DMSO)`，`VeryTightOpt + VeryTightSCF + DefGrid3`；
- 屏蔽生产候选：`revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore/GIAO/CPCM(DMSO)`；
- 稳健性交叉检查：`PBE0/pcSseg-3/GIAO/CPCM(DMSO)`；
- J 耦合候选：`PBE/pcJ-3`，须经独立 J 方法门控后启用。

`RI-MP2`、`TightPNO DLPNO-MP2`、`AUTOCI-CCSD(T)/pcSseg-4` 和 `VPT2` 是默认关闭的受限基准分支。它们更昂贵，但不自动等于更接近溶液实验；只有在资源、数值稳定性、同方法参考和代表集验证全部通过后才能人工启用。详细依据见 [docs/METHOD_RATIONALE.md](docs/METHOD_RATIONALE.md)。

## 已准备与尚未完成

已准备：

- 2 个化合物 × 5 个固定显式溶剂场景 × 20 个起点，共 200 个 GOAT 任务；
- 0、1、2、3 个显式 DMSO，以及 `2 DMSO + 1 H₂O` 痕量水敏感性分支；
- 输入、坐标、原子映射、同位素策略和 SHA-256 provenance；
- 同方法 TMS 参考工作流；
- SLURM/PBS Pro 数组调度、节点本地 scratch、状态与收集工具；
- 苯甲酰胺校准目标和 `(S)`-脯氨酰胺留出检查目标；
- 计算谱与仪器外观层分离的后处理接口。

尚未完成：

- 真实 Linux HPC 上的 ORCA 6.1.1 执行；
- GOAT 构象、连接关系、手性、质子转移和显式溶剂簇的人工批准；
- PBE0 优化/频率、虚频检查和 QRRHO/Boltzmann 审批；
- TMS 参考与分析物 shielding/J 的正常终止和哈希锁定；
- 使用原始实验数据完成的独立留出验证；
- 最终 `¹H`、`¹³C` 预测谱和方法不确定度报告。

## 服务器快速开始

将内层 ASCII 目录 `ORCA_61_NMR_ULTRA_2CMPD` 上传到 Linux 服务器。服务器需自行安装并授权 ORCA 6.1.1；本仓库不包含 ORCA、MPI、许可证、账号口令或私钥。

```bash
cd ORCA_61_NMR_ULTRA_2CMPD
python3 -m pip install --user -r requirements.txt
chmod +x hpc/submit/*.sh hpc/bin/*.py hpc/tests/*.sh

# 必须先改成本站点的 ORCA、调度器和 scratch 设置
cp hpc/config/site_config.json.example hpc/config/site_config.json
$EDITOR hpc/config/site_config.json
# strict_hash=true：下面的 64 位哈希必须写入 orca_sha256
sha256sum /absolute/path/to/orca

# 输入审计与 HPC 预检
python3 scripts/audit_explicit_goat_preparation.py
python3 scripts/build_workflow_profile.py
PROFILE_MANIFEST=data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/GOAT_explicit_job_manifest.csv
python3 hpc/bin/preflight.py \
  --package-root "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --stage goat \
  --scheduler slurm

# 只渲染前两个 SLURM 数组任务；默认不提交
export HPC_SCHEDULER=slurm
hpc/submit/submit_goat.sh \
  "$PROFILE_MANIFEST" \
  --array-range 1-2 \
  --max-concurrent 2

# 两个含水场景也必须各自先做 pilot（仍默认只渲染、不提交）
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 21-22 --max-concurrent 2
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 51-52 --max-concurrent 2
```

人工检查生成的脚本后，才在同一命令末尾添加 `--submit`。PBS Pro 将 `HPC_SCHEDULER` 和预检参数改为 `pbs`。完整服务器步骤见 [docs/SERVER_QUICKSTART.md](docs/SERVER_QUICKSTART.md)，各阶段门控见 [docs/WORKFLOW.md](docs/WORKFLOW.md)。

每一阶段完成后，必须使用 `hpc/bin/collect_results.py --require-all --materialize-job-directories` 逐文件核验 attempt 台账并生成 job-local receipt，之后才能进入下一阶段。VPT2 必须按 `vpt2_force_field → 哈希核验/物化 → vpt2_nmr` 执行；第二阶段不会只凭现有 `.hess`/`.vpt2` 文件放行。

## 一万核服务器的正确用法

不要把 10,000 核分给一个 ORCA 任务。单个任务受 `config/resource_policy.json` 的 rank 上限约束；大规模资源应用于构象、显式溶剂场景、核种、方法和 VPT2 位移的作业数组并发。并发数还必须服从站点内存、scratch、许可证和调度策略。

## 关键目录

```text
01_goat_explicit/       200 个已生成的 ULTRA 母清单输入；默认 practical 从中选取 70 个
config/                 方法、条件、资源和显式溶剂策略
data/                   输入 manifests、场景注册表和审计结果
docs/                   工作流、服务器操作和方法依据
hpc/                    SLURM/PBS Pro 调度、恢复和收集工具
jobs/                   人工审批 gate 模板；默认 fail-closed
reference_jobs/         TMS 同方法参考任务
scripts/                构建器、审计器、解析和谱图工具
validation/             实验峰目标、数据 provenance 与验证政策
deliverables/           只放最终审核通过的表、图和报告
```

## 结果命名规则

最终交付必须同时保留：

- 原始 ORCA 输出及哈希；
- 按原子列出的 shielding 和 J；
- 相同显式溶剂计量内的构象权重；
- `CLEAN_CALCULATED` 峰表/谱；
- 可选 `SYNTHETIC_ACQUISITION_APPEARANCE` 谱；
- 方法分支差异、实验误差指标和未通过项。

`VALIDATED` 只能在所需计算、同方法参考、留出指标和具名复核人齐全后写入。当前状态是 `DRAFT_NOT_VALIDATED`。

---

## English summary

This is a validation-ready, fail-closed ORCA 6.1.1 workflow for `¹H` and `¹³C` NMR of benzamide (`F3_005`) and neutral `(S)`-prolinamide (`F3_015`) in DMSO-d₆. It contains 200 deterministic GOAT inputs, high-level shielding branches, same-method TMS references, provenance hashes, and SLURM/PBS Pro tooling.

The package has **not** been run on the target HPC system and has **not** completed independent experimental validation. Its production candidate is `revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore`; `PBE0/pcSseg-3` is a cross-check. RI-MP2, DLPNO-MP2, AUTOCI-CCSD(T), and VPT2 are gated diagnostics, not automatically superior production methods. Conformers may be Boltzmann averaged only within an identical explicit-solvent composition. Noise, baseline, solvent peaks, intensity, and linewidth belong to a separately labelled synthetic acquisition layer and are not first-principles observables.

### Default practical profile with explicit trace water

`PRACTICAL_HIGH_ACCURACY` is now the declared default daily profile; `ULTRA_DIAGNOSTIC` remains available for adaptive escalation. The practical profile selects three fixed-composition branches for benzamide and four for `(S)`-prolinamide, with ten deterministic GOAT starts per branch (70 GOAT jobs total). For each compound, `DMSO_02_W01_trace_water_sensitivity` (`2 DMSO + 1 H2O`) is mandatory and is compared only with the matched `DMSO_02_W00_primary` dry control. It remains separate from the main prediction, and cross-stoichiometry Boltzmann mixing is forbidden.

The spectrum pipeline preserves a clean calculated layer and can additionally render a separately labelled synthetic acquisition appearance containing empirical noise, baseline, DMSO/water components, peak heights and linewidths. For benzamide `F3_005`, figure-specific single-sample conditioned profiles use the audited `1-1-H/C` appearance descriptors without fitting or moving any calculated analyte chemical shift. See [docs/PRACTICAL_HIGH_ACCURACY.md](docs/PRACTICAL_HIGH_ACCURACY.md).

### 2026-08-10 release-critical scientific controls

- The opt/frequency input evaluates ORCA QRRHO thermochemistry at every nucleus-specific acquisition temperature in `config/nucleus_conditions.csv`, using the same isotope-aware Hessian. `1H` and `13C` weights are stored separately and are never silently reused across temperatures.
- `AUTOCI-CCSD(T)/pcSseg-4` is named and used only as a gas-phase diagnostic. Activation requires a hash-bound ORCA 6.1.1 Linux preflight, the `orca_autoci` helper, and collection of `*.autociresp_*` response outputs.
- Correlated-response diagnostics lock `Z_Tol=1e-12`; the AutoCI branch also locks `STol=1e-10` and `EPRNMR Tol=1e-12`. A branch approval token is insufficient by itself: every named numerical, reference, resource and scope review in `jobs/high_level_manual_opt_in_gate.json` must read `APPROVED_AFTER_REVIEW`.
- VPT2 preserves explicit deuterium masses in both steps. An analyte-only VPT2 result cannot enter the main chemical shift. A reviewed composite requires matched static and vibrational corrections for both analyte and TMS, with all inputs, outputs, receipts and hashes verified.
- The VPT2 force-field and property steps use explicit 0.05 displacements and require `vpt2_ff.hess`, `vpt2_ff.vpt2`, and the property-step `vpt2_nmr.vpt2`. Starting-geometry, soft-mode/resonance, and displacement-sensitivity reviews are mandatory; VPT2 is not presented as a universal improvement.
- Benzamide ortho/meta aromatic H and C pairs are averaged only through the pre-declared identity-locked solution-symmetry layer in `config/solution_equivalence_groups.json`; raw atom-resolved values remain auditable. The two amide NH protons are not symmetry-merged.
- The J benchmark is generated by `scripts/build_j_benchmark_jobs.py`, but remains `enabled=false` until the separate, machine-bound `jobs/j_benchmark_execution_gate.json` is approved. Both PBE/pcJ-3 and PBE0/pcJ-2 lock `Tol=1e-10`, `Z_Tol=1e-12`, and the 6.05 Å pair screen. Routine carbon output remains proton-decoupled `13C{1H}`.
- Validation is controlled by pre-locked signal IDs and atom/symmetry provenance. Hungarian minimum-error matching is auxiliary only. The strongest possible status is `VALIDATED_FOR_THIS_TWO_COMPOUND_PROTOCOL_AND_LOCKED_CONDITIONS`, never a general claim for drug-like molecules.
