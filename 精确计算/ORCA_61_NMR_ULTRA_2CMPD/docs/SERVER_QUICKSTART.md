# Linux 服务器快速运行手册 / Server quick start

## 1. 前置条件

服务器必须具备：

- 64 位 Linux；
- 已授权的 ORCA 6.1.1 Linux 发行包及其匹配的 MPI；
- ORCA 所需 xTB 组件或站点确认的兼容配置；
- Python 3.10 或更高版本；
- SLURM 或 PBS Pro；
- 充足的节点本地 scratch；
- 允许的 ORCA 并发数和许可证策略。

本包不包含 ORCA、MPI、xTB 许可证、账号、密码、SSH key 或其他秘密。不要把这些内容提交到 Git。

## 2. 上传与目录名

推荐直接从 GitHub 克隆，然后进入内层 ASCII 目录，避免部分旧集群对中文路径处理不一致：

```bash
git clone git@github.com:zqk-7k/DK_Computation.git
cd DK_Computation/精确计算/ORCA_61_NMR_ULTRA_2CMPD
```

若集群不能稳定处理中文路径，可在本地自行把内层 `ORCA_61_NMR_ULTRA_2CMPD` 目录打成 tar.gz 后上传；仓库本身不提交预制 tar 包。解包后直接进入该 ASCII 目录。

不要在 Windows 上生成 Linux ORCA 输出后再把它们当作服务器结果。生产计算必须使用服务器上的 Linux ORCA 6.1.1。

## 3. 安装 Python 依赖和设置权限

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt

chmod +x hpc/submit/*.sh hpc/bin/*.py hpc/tests/*.sh
```

`.venv/`、ORCA 二进制和计算 scratch 不应进入版本控制。

## 4. 配置站点

```bash
cp hpc/config/site_config.json.example hpc/config/site_config.json
$EDITOR hpc/config/site_config.json
```

至少检查：

- `orca_executable`：ORCA 主程序绝对路径；
- `xtb_executable`：与 ORCA 发行包匹配的 xTB 路径；
- `module_commands`：站点 module 加载命令；
- SLURM `account/partition/qos` 或 PBS `project/queue`；
- `scratch_base`；
- `orca_sha256`：Linux ORCA 主程序的准确 64 位十六进制 SHA-256；`strict_hash=true` 时为必填；
- AutoCI 若获人工批准，填写 `helper_executables.orca_autoci`。

在服务器上生成哈希并复制到配置：

```bash
sha256sum /absolute/path/to/orca
```

ORCA 主程序用绝对路径直接启动：

```text
/absolute/path/to/orca input.inp
```

不要在外层再套 `mpirun` 或 `mpiexec`。ORCA 根据 `%pal` 自行启动 MPI 子进程。

## 5. 输入审计

```bash
python3 scripts/audit_explicit_goat_preparation.py
```

该审计预期检查 200 个 ULTRA 母清单输入、坐标、原子/同位素映射及 SHA-256。任何失败都应先修复，不要通过删除 manifest 行或关闭审计绕过。随后生成/刷新默认 practical 子清单：

```bash
python3 scripts/build_workflow_profile.py
PROFILE_MANIFEST=data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/GOAT_explicit_job_manifest.csv
PROFILE_SCENARIOS=data/workflow_profiles/PRACTICAL_HIGH_ACCURACY/explicit_scenario_registry.csv
```

默认 practical 包含 70 条；200 条母清单只用于 `ULTRA_DIAGNOSTIC`。

## 6. 调度器预检

SLURM：

```bash
python3 hpc/bin/preflight.py \
  --package-root "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --stage goat \
  --scheduler slurm
```

PBS Pro：

```bash
python3 hpc/bin/preflight.py \
  --package-root "$PWD" \
  --manifest "$PROFILE_MANIFEST" \
  --stage goat \
  --scheduler pbs
```

预检必须确认 ORCA/xTB 路径、调度命令、相对路径、输入哈希、`%pal`、GOAT `NWorkers`、`%maxcore` 和申请资源一致。

## 7. 两任务 pilot

SLURM：

```bash
export HPC_SCHEDULER=slurm
hpc/submit/submit_goat.sh \
  "$PROFILE_MANIFEST" \
  --array-range 1-2 \
  --max-concurrent 2
```

PBS Pro：

```bash
export HPC_SCHEDULER=pbs
hpc/submit/submit_goat.sh \
  "$PROFILE_MANIFEST" \
  --array-range 1-2 \
  --max-concurrent 2
```

默认行为只生成调度脚本，不会提交。打开 `state/generated_submissions/` 中的新脚本，检查账户、队列、CPU、内存、walltime、数组范围、scratch 和 ORCA 路径。确认后在相同命令末尾添加：

```text
--submit
```

普通 pilot 通过后，还必须分别渲染并提交两个含水场景的 pilot：

```bash
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 21-22 --max-concurrent 2
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 51-52 --max-concurrent 2
```

检查调度脚本后再追加 `--submit`。这些范围对应两个化合物各自的 `2 DMSO + 1 H2O` 分支。

## 8. 扩大批次

所有 pilot 成功后，先扩大普通分支并分别跑完整的 10-task 含水分支：

```bash
hpc/submit/submit_goat.sh \
  "$PROFILE_MANIFEST" \
  --array-range 1-10 \
  --max-concurrent 4 \
  --submit
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 21-30 --max-concurrent 4 --submit
hpc/submit/submit_goat.sh "$PROFILE_MANIFEST" --array-range 51-60 --max-concurrent 4 --submit
```

稳定后删除 `--array-range` 才会覆盖默认 practical 的全部 70 条；已经成功的任务会按已核验状态跳过：

```bash
hpc/submit/submit_goat.sh \
  "$PROFILE_MANIFEST" \
  --max-concurrent 16 \
  --submit
```

GOAT 每任务申请 8 CPU 和 24 GiB。初次全量建议并发 16；请根据每节点内存、scratch I/O、许可证并发和队列政策调整。只有 pilot 和扩展批次都稳定时才逐步升到 32 或 60。60 并发需要 480 CPU，并应至少有约 1.44 TiB 可调度内存。提高 `--max-concurrent` 只增加并行任务数，不提高单个计算的理论精度。

## 9. 状态、失败与收集

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
```

结果位于 `results/hpc/<stage>/<job_id>/attempt_NNN/`，状态位于 `state/tasks/<stage>/`。`--materialize-job-directories` 会重新核验源台账，把无冲突结果复制回每个 job 目录，并写入跨阶段 receipt。VPT2 必须先完成并物化 `vpt2_force_field`，之后 `vpt2_nmr` 才能通过依赖预检。重提同一 manifest 会跳过已重新核验成功的任务，只处理失败或未完成任务。

若发现 stale lock，必须先用 `squeue` 或 `qstat` 确认对应调度任务已经结束，再使用 `hpc/bin/unlock_stale_task.py`。不要在作业仍运行时删除锁。

## 10. 下一阶段不能自动串联

GOAT 完成不等于可以直接计算 NMR。必须先导出兼容状态、生成草稿 gate，并由具名人员审查：

```bash
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

之后按 [WORKFLOW.md](WORKFLOW.md) 完成人工批准、PBE0-D4 opt/freq、TMS、shielding、J、构象平均和验证。

## 11. 一万核资源策略

将大规模资源用于作业数组，而不是单个超大 MPI 作业：

- GOAT 起点可并发；
- 不同固定溶剂计量可并发，但结果不能混合平均；
- opt/freq 构象可并发；
- revDSD 与 PBE0 分支可并发；
- 不同核种和 VPT2 位移可并发；
- 单任务 rank 不超过 `config/resource_policy.json` 的绝对上限。

10,000 个可用核不表示应同时占满。先以 pilot 测得实际内存、scratch 和 walltime，再计算安全数组并发。

## 12. 结束检查

在宣布任何结果前确认：

- 所有相关任务正常终止；
- 输入与输出 hashes 可追溯；
- 参考和分析物 method fingerprint 完全一致；
- 无显著虚频、错误物种、质子转移或手性翻转；
- 只在相同显式计量内平均；
- clean 计算谱和 synthetic 外观谱分开；
- 留出验证和具名复核已完成。

当前仓库仅满足“输入已准备”，并不满足上述最终检查。

---

## English quick runbook

Upload the inner ASCII directory to Linux, install the Python requirements, configure `hpc/config/site_config.json`, and run the preparation audit, practical-profile materialization, and scheduler preflight. Render dry and trace-water pilot arrays without `--submit`, inspect the generated scheduler scripts, and then submit them. Expand each wet branch to ten tasks and finally run the 70-job practical manifest only after every pilot succeeds. The 200-job master manifest is reserved for `ULTRA_DIAGNOSTIC`.

Use a large machine for array concurrency, not as a single 10,000-rank ORCA job. Do not chain GOAT directly into NMR properties: connectivity, stereochemistry, cluster integrity, stationary points, thermochemistry, same-method TMS references, and every manual gate must be reviewed first. This repository currently contains prepared inputs, not completed HPC results or validated spectra.
