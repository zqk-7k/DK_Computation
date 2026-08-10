# 工作流与人工门控 / Workflow and manual gates

## 1. 总体流程

```text
结构与条件锁定
  → 固定溶剂计量的 GOAT 构象搜索
  → 人工构象/连接/手性门控
  → PBE0-D4 优化与频率
  → 虚频、簇完整性和 QRRHO/Boltzmann 门控
  → 同方法 TMS 参考
  → revDSD 生产候选 + PBE0 交叉检查 shielding
  → 受控 J 耦合与可选高阶基准
  → 同计量构象平均和参考定标
  → clean calculated spectrum
  → 独立 synthetic acquisition appearance
  → 苯甲酰胺校准 + 脯氨酰胺留出验证
```

所有构建器只生成输入，不启动 ORCA。所有人工 gate 默认 `DRAFT_NOT_APPROVED`；不得用自动脚本把草稿批量改成批准状态。

## 2. 阶段 0：身份与实验条件锁定

使用以下文件：

- `registry/Figure3_52_ORCA_reference_state_registry.csv`
- `config/compound_queue.csv`
- `config/Figure3_52_proposed_NMR_conditions.csv`
- `config/nucleus_conditions.csv`

必须确认分子身份、净电荷、多重度、立体化学、质子化态、盐型、溶剂、核种、温度和参考标度。`F3_015` 当前假设为中性 `(S)`-脯氨酰胺游离碱；它的实验元数据来自处理后图像，未达到原始 FID 级完整性。任何身份或条件变化都要求重新生成输入和哈希。

## 3. 阶段 1：固定计量 GOAT 搜索

场景由 `config/microsolvation_policy.json` 锁定：

| 场景 | 显式环境 | 用途 |
|---|---|---|
| `DMSO_00_W00` | 无显式溶剂 | 连续介质对照 |
| `DMSO_01_W00` | 1 DMSO | 第一溶剂壳敏感性 |
| `DMSO_02_W00` | 2 DMSO | 苯甲酰胺主要候选 |
| `DMSO_03_W00` | 3 DMSO | 脯氨酰胺扩展候选 |
| `DMSO_02_W01` | 2 DMSO + 1 H₂O | 痕量水敏感性，不代表样品水活度 |

ULTRA 母清单中每个化合物、每个场景有 20 个确定性起点，共 200 个任务。重建母输入时运行：

```bash
python3 scripts/prepare_explicit_goat_jobs.py \
  --starts-per-scenario 20 \
  --force-inputs
python3 scripts/audit_explicit_goat_preparation.py
```

日常生产默认从母清单物化 `PRACTICAL_HIGH_ACCURACY`：苯甲酰胺 3 个场景、脯氨酰胺 4 个场景，每场景 10 个起点，共 70 个任务。200 个任务只在自适应门控升级到 `ULTRA_DIAGNOSTIC` 后使用。

```bash
python3 scripts/build_workflow_profile.py
PROFILE_MANIFEST=data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/GOAT_explicit_job_manifest.csv
PROFILE_SCENARIOS=data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/explicit_scenario_registry.csv
```

先运行普通 2-task pilot，并对清单行 21–22 与 51–52 的两个含水分支分别做 pilot；再扩展到每个含水分支 10 个任务，确认稳定后才提交 practical 全数组。

GOAT 完成后，使用 HPC 收集与桥接工具：

```bash
python3 hpc/bin/status.py \
  --package-root "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --stage goat

python3 hpc/bin/collect_results.py \
  --package-root "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --stage goat \
  --require-all \
  --materialize-job-directories

python3 hpc/bin/export_workflow_states.py \
  --package-root "$PWD" \
  --stage goat \
  --manifest "$PROFILE_MANIFEST" \
  --output-state-root logs/goat_explicit_batch/jobs \
  --compat-manifest data/GOAT_manifest_gate_compatible.csv

python3 scripts/draft_goat_selection_gate.py \
  --project "$PWD" \
  --manifest data/GOAT_manifest_gate_compatible.csv \
  --scenario-registry "$PROFILE_SCENARIOS" \
  --state-root logs/goat_explicit_batch/jobs
```

人工审查每个保留构象的连接关系、原子顺序、净电荷、多重度、质子转移、显式溶剂完整性及 `F3_015` 的 S 手性。只有通过后，复制/改名为 `goat_ensemble_selection_approved.json` 并按模板填写 `APPROVED_FOR_OPT_FREQ`、复核人、时间和理由。

## 4. 阶段 2：PBE0-D4 优化与频率

```bash
python3 scripts/build_dft_jobs.py \
  --project "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --scenario-registry "$PROFILE_SCENARIOS"

python3 hpc/bin/normalize_manifest_paths.py \
  --package-root "$PWD" \
  --input data/DFT_explicit_job_manifest.csv \
  --output data/DFT_explicit_job_manifest_portable.csv \
  --stage optfreq

hpc/submit/submit_optfreq.sh \
  data/DFT_explicit_job_manifest_portable.csv \
  --max-concurrent 8
```

上面的提交命令默认只是 dry-run；审核后追加 `--submit`。`--max-concurrent` 是站点并发控制，不改变单任务理论层级。

完成后桥接状态并生成 minimum/Boltzmann 草稿：

```bash
python3 hpc/bin/collect_results.py \
  --package-root "$PWD" \
  --manifest data/DFT_explicit_job_manifest_portable.csv \
  --stage optfreq \
  --require-all \
  --materialize-job-directories

python3 hpc/bin/export_workflow_states.py \
  --package-root "$PWD" \
  --stage optfreq \
  --manifest data/DFT_explicit_job_manifest_portable.csv \
  --output-state-root logs/opt_freq_serial/jobs \
  --compat-manifest data/DFT_manifest_gate_compatible.csv

python3 scripts/draft_minimum_boltzmann_gate.py \
  --project "$PWD" \
  --manifest data/DFT_manifest_gate_compatible.csv \
  --state-root logs/opt_freq_serial/jobs
```

批准前必须核对：正常终止、无显著虚频、几何仍是目标物种、显式簇未解离、`F3_015` 手性未翻转、热化学可解析、温度正确。审批状态为 `APPROVED_FOR_NMR`。

## 5. 阶段 3：QRRHO 与构象平均

权重必须按以下键分别归一化：

```text
figure_id + scenario_id + nucleus + acquisition_temperature
```

禁止把 0、1、2、3 DMSO 或 `2 DMSO + 1 H₂O` 用普通 Boltzmann 权重混在一起。这些模型含有不同分子数和不同标准态；除非另行建立并验证溶剂化学势/活度模型，否则只能作为并列的模型分支比较。

`¹H` 和 `¹³C` 若采集温度不同，必须分别计算权重。若剔除构象，必须仅用剩余构象重新归一化。

## 6. 阶段 4：TMS 同方法参考

先准备参考输入：

```bash
python3 scripts/prepare_ultra_reference_jobs.py
```

运行顺序为：

1. `TMS_DMSO__opt_freq`；
2. 检查正常终止、零虚频和 `tms_opt.xyz`；
3. 锁定坐标 SHA-256；
4. 仅启用与分析物完全相同的 shielding 方法分支；
5. 运行 reference shielding；
6. 将输出、方法 fingerprint、坐标和输出哈希写入 `config/reference_job_registry.json`；
7. 具名复核后改为批准状态。

第 2 步人工审查通过后，使用 `scripts/activate_tms_reference_jobs.py`。它会生成只包含本次批准 shielding 分支的 `data/reference_shielding_active_manifest.csv`；不要把同时含 opt/freq 与 shielding 的总表直接当作激活清单。每个参考阶段完成后同样必须使用 `collect_results.py --materialize-job-directories --require-all`。

目标与参考必须在方法、基组、溶剂模型、冻结核设置以及 RI/局域相关近似上完全一致。不得用 PBE0-TMS 标定 revDSD 分析物，也不得用经验常数替代缺失的同方法参考却仍宣称 ab initio 绝对定标。

## 7. 阶段 5：shielding 与 J

在 minimum gates 和参考 registry 均批准后生成属性输入：

```bash
python3 scripts/build_property_jobs.py \
  --project "$PWD" \
  --workflow-profiles config/workflow_profiles.json \
  --workflow-profile PRACTICAL_HIGH_ACCURACY
```

该构建器同时生成：

- `production_revdsd_pcsseg3`：生产候选；
- `crosscheck_pbe0_pcsseg3`：全构象交叉检查；
- J 输入：只有 `jobs/j_method_benchmark_gate.json` 通过后才能生成/执行。

对生成的 manifest 使用 `normalize_manifest_paths.py`，随后依次使用 `submit_reference_shielding.sh`、`submit_shielding.sh` 和 `submit_jcoupling.sh`。不要把不同 method fingerprint 的 shielding 混入同一次构象平均。

可选高阶分支通过：

```bash
python3 scripts/build_high_level_benchmarks.py \
  --project "$PWD" \
  --gate jobs/F3_005/minimum_boltzmann_audit/minimum_boltzmann_approved.json \
  --gate jobs/F3_015/minimum_boltzmann_audit/minimum_boltzmann_approved.json
```

其任务仍以 `enabled=false` 和 `readiness=manual_opt_in` 写出。只有 `jobs/high_level_manual_opt_in_gate.json` 完成资源、数值稳定性、同方法参考和用途审查后才能逐分支启用。

VPT2 是两个独立调度阶段：先提交 `submit_vpt2_force_field.sh`，成功后用 `collect_results.py --stage vpt2_force_field --materialize-job-directories --require-all` 生成哈希绑定 receipt，再提交 `submit_vpt2_nmr.sh`。缺失、冲突或被改动的 `.hess`/`.vpt2` 及其源台账都会阻止第二阶段。J 方法代表集由 `build_j_benchmark_jobs.py` 生成，并使用 `submit_j_method_benchmark.sh`；其结果只用于方法选择。

## 8. 阶段 6：化学位移与谱图

化学位移的基本关系为：

```text
δ_i = σ_reference − σ_i
```

每个 method fingerprint 使用自己的同方法 TMS。构象平均可对同场景内的 shielding 先平均，再换算位移；无论采用何种顺序，都应保留每个原子、每个构象、每个权重和参考值的表格。

输出必须分为：

- `RAW_ORCA_PROPERTIES`：未经经验修正的 shielding/J；
- `CLEAN_CALCULATED`：构象平均、参考定标和自旋模拟后的计算谱；
- `SYNTHETIC_ACQUISITION_APPEARANCE`：加入溶剂、噪声、基线和线宽模型的外观谱。

外观层只能改变显示和采集模型，不能借实验峰逐峰移动计算峰。交换性 NH 必须单独报告误差；默认不把其计算 J 强制表现为相干精细结构。

## 9. 阶段 7：验证与晋级

`F3_005` 是方法校准化合物，`F3_015` 是留出效果检查。目标和限制在：

- `validation/experimental_peak_targets.csv`
- `validation/validation_policy.json`
- `validation/benzamide_raw_fid_provenance.json`

主要目标是非交换 `¹H` MAE ≤ 0.20 ppm、`¹³C` MAE ≤ 2.5 ppm，同时人工检查峰数、顺序和多重性。它们是晋级目标，不是当前达到的结果，也不是对新分子的保证。`F3_015` 目前缺少完整原始采集元数据，所以即使数值通过也只能暂定通过。

严禁：

- 用 `F3_015` 的实验峰逐峰调参后再称其为留出集；
- 复制实验噪声轨迹；
- 混合不同显式溶剂计量；
- 将 synthetic 外观称为实验谱；
- 只报告最好方法而隐藏失败分支。

## 10. 审计与恢复

每个阶段都应保留输入、坐标、输出、状态 JSON、文件清单和 SHA-256。失败任务进入新的 attempt，不覆盖旧 attempt。只有返回码为 0、输出包含 `ORCA TERMINATED NORMALLY`、必要产物存在且输入哈希匹配，任务才算完成。

不要把 ORCA 二进制、许可证、服务器账号、Git 凭据或私钥写入仓库。

---

## English workflow summary

The workflow is deliberately staged and fail-closed. Lock identity and acquisition conditions; run GOAT separately for each fixed explicit-solvent composition; manually approve connectivity, protonation, stereochemistry, and cluster integrity; optimize and compute frequencies at PBE0-D4/def2-TZVPP; approve stationary points and nucleus-specific QRRHO weights; compute a same-method TMS reference; run the revDSD production candidate and PBE0 cross-check; then generate clean calculated spectra and a separately labelled synthetic acquisition appearance.

Ordinary Boltzmann averaging across different explicit-solvent stoichiometries is forbidden. The same isotope-aware Hessian is evaluated at every locked nucleus acquisition temperature, and the resulting nucleus-specific QRRHO weights are authoritative. RI-MP2, DLPNO-MP2, gas-phase AUTOCI-CCSD(T), and VPT2 remain machine-bound diagnostics until their manual gates are reviewed. Every named `scientific_reviews` item must be `APPROVED_AFTER_REVIEW`; a parent or branch token alone cannot enable execution. AutoCI locks `STol=1e-10`, response `Tol=1e-12`, and `Z_Tol=1e-12`. VPT2 requires the exact force-field pair `vpt2_ff.hess`/`vpt2_ff.vpt2`, the property-step `vpt2_nmr.vpt2`, and reviewed starting-geometry, resonance, soft-mode, and displacement sensitivity. Before approving the production J method, build the paired representative inputs with `python3 scripts/build_j_benchmark_jobs.py --gate <F3_005-approved-gate> --gate <F3_015-approved-gate>`. Input generation does not authorize execution: `jobs/j_benchmark_execution_gate.json` must bind both minimum-gate hashes, the resource policy, the target ORCA 6.1.1 preflight and a named reviewer. After execution, review method and 6.05 Å screen sensitivity in `jobs/j_method_benchmark_gate.json`. VPT2 cannot enter the main delta without matched analyte/TMS static and vibrational corrections. Validation uses fixed signal IDs and atom/symmetry provenance; Hungarian matching is auxiliary only. The only promotion token is `VALIDATED_FOR_THIS_TWO_COMPOUND_PROTOCOL_AND_LOCKED_CONDITIONS`, and it requires complete calculations, hashes, references, held-out metadata/metrics, and a named reviewer.
