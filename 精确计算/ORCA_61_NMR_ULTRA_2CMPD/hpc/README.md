# 两化合物 ORCA 6.1.1 HPC 调度层

本目录仅负责在 Linux SLURM/PBS Pro 上安全地运行 manifest。完整流程见
[`../docs/SERVER_QUICKSTART.md`](../docs/SERVER_QUICKSTART.md)。当前准备好的第一阶段是
`data/GOAT_explicit_job_manifest.csv`：2 个化合物、10 个固定显式溶剂情景、每个情景 20 个独立起点，共 200 个任务。

所有任务均为 **CALCULATED MODEL / NOT EXPERIMENTAL DATA**。当前仓库只包含输入，未包含生产 ORCA 结果。

## 站点配置

复制并编辑 `config/site_config.json.example` 为本机专用且被 Git 忽略的 `config/site_config.json`。必须填写 Linux ORCA 6.1.1 可执行文件的绝对路径及其准确 SHA-256、调度账户/分区和可选 scratch。`strict_hash=true` 时 SHA-256 缺失或不匹配都会阻止提交。ORCA 主程序必须用绝对路径直接启动，不能在外面再套 `mpirun`。AutoCI 分支还必须在 `helper_executables` 中绑定 `orca_autoci`，或把它安装在 ORCA 主程序旁。

资源由 `config/stage_resources.json` 控制。生产输入使用：

- GOAT：8 ranks，`%maxcore 2000`；
- opt/freq：16 ranks，`%maxcore 6000`；
- shielding：16 ranks，`%maxcore 8000`；
- J coupling：16 ranks，`%maxcore 6000`；
- TMS shielding 预留到 16 ranks × 12 GB/rank。

一万核服务器应通过数组并发利用，而不是把一万个核分给一个 ORCA 作业。

## 预检与提交

```bash
python3 hpc/bin/preflight.py \
  --package-root "$PWD" \
  --manifest data/GOAT_explicit_job_manifest.csv \
  --stage goat \
  --scheduler slurm

# 首先只渲染，不提交
export HPC_SCHEDULER=slurm
hpc/submit/submit_goat.sh \
  data/GOAT_explicit_job_manifest.csv \
  --array-range 1-2 --max-concurrent 2

# 人工检查生成的脚本后，才允许真实提交
hpc/submit/submit_goat.sh \
  data/GOAT_explicit_job_manifest.csv \
  --max-concurrent 100 --submit
```

PBS Pro 把 `HPC_SCHEDULER` 和 preflight 的 scheduler 改为 `pbs`。

每个 array task 只运行一个 ORCA 输入。runner 会检查输入哈希、CPU/内存配额、正常终止标志和阶段必需产物，并把每次尝试写入 `results/hpc/` 与 `state/tasks/`。失败任务不会覆盖旧尝试。

每个阶段完成后必须运行 `collect_results.py --materialize-job-directories --require-all`。该命令逐文件核验 attempt 台账，然后无冲突地把结果物化回 job 目录并写入依赖 receipt；后续阶段不得只相信文件存在。

## 必须保持的门控

GOAT 完成后不能直接串联所有后续阶段。必须依次：收集结果、导出兼容状态、人工批准构象选择、DFT opt/freq、人工批准无显著虚频的最低点、同方法 TMS、shielding、J、解析与验证。不同 DMSO/H2O 显式计量必须始终作为独立模型分支，不能做普通 Boltzmann 混合。

RI-MP2、DLPNO-MP2、AUTOCI-CCSD(T) 和 VPT2 清单默认 `enabled=false`，只能在人工资源/科学门控后启用。
