# ORCA 52 化合物 Linux 超算调度层

本目录把便携计算包转换为可审计的 SLURM 或 PBS Pro 作业数组。当前已可直接提交的生产阶段是 **GOAT 构象搜索**，其 1,120 条任务清单为：

`data/GOAT_explicit_job_manifest.csv`

`optfreq`、`reference_optfreq`、`reference_shielding`、`shielding` 和 `jcoupling` 已有相同的调度模板与资源策略；这些阶段必须先通过前一阶段的人工科学门控，并由工作流脚本生成对应的便携清单后再提交。调度层不会跳过结构、构象、虚频、Boltzmann 权重或参考物审批。

> 所有输出均为 **CALCULATED MODEL / NOT EXPERIMENTAL DATA**。ORCA 不能生成实验 HRMS 碎裂图；本包中的质谱部分是分子式、精确质量与同位素包络预测，不应称为实测质谱。

## 1. 上传后只需修改一个站点文件

在超算登录节点解压后进入包根目录，然后编辑：

`hpc/config/site_config.json`

至少设置：

- `orca_executable`：ORCA 6.1.1 Linux 可执行文件的**绝对路径**；
- `xtb_executable`：GOAT 所用、与该 ORCA 发行包匹配的 xTB 路径；
- `module_commands`：集群需要的 module 命令，可为空；
- SLURM 的 `account/partition/qos`，或 PBS 的 `project/queue`；
- 可选的 `scratch_base`。留空时优先使用 `SLURM_TMPDIR`，其次使用 `TMPDIR`。

示例在 `hpc/config/site_config.json.example`。不要把 Windows 版 ORCA 上传到 Linux 节点运行。

ORCA 并行作业必须用 ORCA 可执行文件的完整路径直接启动。这里的任务运行器执行：

```text
/absolute/path/to/orca input.inp
```

它不会调用 `mpirun` 或 `mpiexec`；ORCA 自己负责启动 MPI 子进程。

## 2. 初次权限和预检

```bash
cd /path/to/ORCA_52_HPC_portable_package
chmod +x hpc/submit/*.sh hpc/bin/*.py hpc/tests/*.sh

python3 hpc/bin/preflight.py \
  --package-root "$PWD" \
  --manifest data/GOAT_explicit_job_manifest.csv \
  --stage goat \
  --scheduler slurm
```

PBS Pro 把最后一项改成 `--scheduler pbs`。预检会检查：Linux、调度器命令、ORCA/xTB 路径、ORCA 二进制哈希（若填写）、1,120 个相对输入路径、输入哈希、`%pal`/`NWorkers`、`%maxcore` 与申请资源是否匹配。报告保存在 `state/preflight/`。

## 3. 先 dry-run 渲染，不提交

SLURM：

```bash
export HPC_SCHEDULER=slurm
hpc/submit/submit_goat.sh data/GOAT_explicit_job_manifest.csv --array-range 1-2 --max-concurrent 2
```

PBS Pro：

```bash
export HPC_SCHEDULER=pbs
hpc/submit/submit_goat.sh data/GOAT_explicit_job_manifest.csv --array-range 1-2 --max-concurrent 2
```

默认只在 `state/generated_submissions/` 生成并记录调度脚本，绝不调用 `sbatch/qsub`。上面的 `1-2` 是安全 pilot，只运行 manifest 的前两行；通过后可先扩到 `--array-range 1-10`，最终正式批次删除 `--array-range` 才会覆盖 1,120 行。范围严格按一基索引检查，不能越过 manifest。请人工检查生成脚本中的账户、分区/队列、CPU、内存、walltime、数组范围和 ORCA 路径。

GOAT 默认是 4 CPU、8 GiB、72 h、最多并发 16 个任务。简单分子在本机也可能需要数小时，大显式簇可能明显超过 24 h；72 h 可减少达到时限后反复重算。若站点单任务上限不同，请在 `stage_resources.json` 中按站点规则调整 walltime。`--max-concurrent` 只改变同时占用的节点资源，不改变计算方法或单任务质量。若集群限制更严，可降到 4–8；不要通过删构象、减起始构型或改理论水平来加速。

## 4. 正式提交

确认 dry-run 脚本后：

```bash
hpc/submit/submit_goat.sh \
  data/GOAT_explicit_job_manifest.csv \
  --max-concurrent 12 \
  --submit
```

`--submit` 会再次运行 preflight；只有通过后才调用 `sbatch` 或 `qsub`。提交命令、脚本哈希、manifest 哈希和调度器返回值保存在同名 `.json` 文件。

每个数组元素只运行一个 ORCA 输入。调度器分配 4 个 CPU 时：

- ORCA 输入中的 `%pal nprocs 4` 与 GOAT `NWorkers 4` 不会超过配额；
- `OMP_NUM_THREADS`、MKL、OpenBLAS、BLIS 和 NumExpr 均强制为 1；
- 不会出现“数组并发 × ORCA 线程 × BLAS 线程”的隐式超额占用。

## 5. 节点本地 scratch、失败和恢复

每个任务先把哈希锁定的输入目录复制到节点本地 scratch，在 scratch 中运行，随后把必需结果原子性复制到：

`results/hpc/<stage>/<job_id>/attempt_NNN/`

状态文件位于：

`state/tasks/<stage>/<job_id>.json`

成功必须同时满足：

1. 进程返回码为 0；
2. 输出含 `ORCA TERMINATED NORMALLY`；
3. 当前阶段要求的文件存在且非空；
4. 输入哈希仍与 manifest 相同。

失败结果不会覆盖旧 attempt。默认保留失败 scratch，具体路径写入状态文件。重新提交相同 manifest 即可恢复：已经完成且重新核验通过的任务会直接标为 `RESUMED_ALREADY_COMPLETE`，只重新计算未完成或失败的任务。

若节点被强制终止，可能遗留锁。先用 `squeue`/`qstat` 确认对应任务已不再运行，再执行：

```bash
python3 hpc/bin/unlock_stale_task.py \
  --package-root "$PWD" \
  --stage goat \
  --job-id F3_001_DMSO_01_primary_s01 \
  --confirm-scheduler-job-is-not-running \
  --minimum-age-seconds 3600
```

解锁操作会写审计记录。不要在任务仍运行时删除锁。

## 6. 查看状态与收集

```bash
python3 hpc/bin/status.py \
  --package-root "$PWD" \
  --manifest data/GOAT_explicit_job_manifest.csv \
  --stage goat

python3 hpc/bin/collect_results.py \
  --package-root "$PWD" \
  --manifest data/GOAT_explicit_job_manifest.csv \
  --stage goat \
  --require-all
```

状态表写入 `state/status/`；完成结果的 SHA256 索引写入 `results/hpc/goat/`。`--require-all` 在任何启用任务未完成时返回非零，适合用作下一阶段的门控。

HPC runner 的状态格式与原 Windows 串行 runner 的旧门控格式不同。GOAT 全部完成后，必须先运行经过哈希验证的桥接器；它不会批准任何人工门控：

```bash
python3 hpc/bin/export_workflow_states.py \
  --package-root "$PWD" \
  --stage goat \
  --manifest data/GOAT_explicit_job_manifest.csv \
  --output-state-root logs/goat_explicit_batch/jobs \
  --compat-manifest data/GOAT_manifest_gate_compatible.csv

python3 scripts/draft_goat_selection_gate.py \
  --project "$PWD" \
  --manifest data/GOAT_manifest_gate_compatible.csv \
  --scenario-registry data/explicit_scenario_registry.csv \
  --state-root logs/goat_explicit_batch/jobs
```

桥接器只接受 `COMPLETED/RESUMED_ALREADY_COMPLETE` 状态，逐文件核验 `FILE_MANIFEST_SHA256.csv`、源输入和 bundle 哈希，再输出旧 validator 可读取的 `status=completed`/`attempts[]` 结构和一份兼容 manifest。生成的 selection gate 仍是 `DRAFT_NOT_APPROVED`，必须人工审查后才能改成批准状态。

## 7. 后续阶段的固定顺序

计算链不能用单个依赖参数盲目串到底。推荐顺序是：

1. `goat` 数组完成并收集；
2. 人工审查构象、连接关系、质子转移、显式溶剂/反离子完整性，并锁定选中 XYZ 哈希；
3. 用工作流脚本生成 opt/freq manifest，再用下述规范化工具把构建器写出的绝对路径转成包内相对路径并重新计算哈希，然后提交 `optfreq`；
4. 审查正常终止、虚频、簇完整性和 QRRHO/Boltzmann 权重；
5. 先计算并批准 `reference_optfreq` 和 `reference_shielding`；
6. 生成并提交分析物 `shielding`；
7. 代表集验证通过后生成并提交 `jcoupling`；
8. 解析 shielding/J，进行自旋模拟和带明确 synthetic 标签的采集外观后处理。

opt/freq 完成后，在建立 minimum/Boltzmann draft gate 前同样桥接：

```bash
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

输入构建器使用 `config/resource_policy.json`：每个 ORCA 输入仍固定为最多 6 个 MPI rank（特殊大簇 4 个），但诚实声明为“每个调度数组元素一个 ORCA 计算”，不虚假宣称整个超算队列全局只有一个化合物运行。全局并发由数组 `%N` 控制。

通用提交命令：

```bash
hpc/submit/submit_optfreq.sh MANIFEST.csv --max-concurrent 6
hpc/submit/submit_reference_optfreq.sh MANIFEST.csv --max-concurrent 2
hpc/submit/submit_reference_shielding.sh MANIFEST.csv --max-concurrent 2
hpc/submit/submit_shielding.sh MANIFEST.csv --max-concurrent 6
hpc/submit/submit_jcoupling.sh MANIFEST.csv --max-concurrent 4
```

下游构建器可能在 `job_directory` 中写入当前登录节点的绝对 POSIX 路径。每次提交新阶段前先生成一个新的便携 manifest（不要原地覆盖源 manifest）：

```bash
python3 hpc/bin/normalize_manifest_paths.py \
  --package-root "$PWD" \
  --input path/to/generated_optfreq_manifest.csv \
  --output data/optfreq_manifest_portable.csv \
  --stage optfreq
```

`--stage` 可取 `optfreq`、`reference_optfreq`、`reference_shielding`、`shielding` 或 `jcoupling`。工具只接受包根目录内的 job/input/coordinate 路径，拒绝无法安全映射的包外路径，补齐相对路径字段并重新计算 input/coordinate（以及适用时的 GOAT bundle）SHA256。随后对生成的 `*_portable.csv` 运行 preflight，再交给对应的 `submit_*.sh`。

第一次仍是 dry-run；人工检查后追加 `--submit`。若前一批作业已经经过人工门控，才可以用 `--dependency JOB_ID` 生成 `afterok` 依赖。

## 8. 下游 manifest 最低字段

CSV 每一行代表一个独立 ORCA 任务，最低需要：

```text
job_id,stage,portable_job_directory,input_file,output_file,expected_input_sha256,enabled,readiness
```

路径必须相对于包根目录，不能包含 `C:\...`、`D:\...` 或 `E:\...`。建议同时记录 `expected_coordinate_sha256`、`required_artifacts`、结构/构象哈希、charge、multiplicity、solvent、temperature、model fingerprint 和审批状态。严格模式下，没有 expected input hash 的任务会被拒绝。

## 9. 离线测试

测试不会启动 ORCA；它使用伪 ORCA 验证 runner、resume、哈希拒绝、CPU 超额拒绝、SLURM/PBS 渲染和“模板中无 mpirun”：

```bash
bash hpc/tests/run_tests.sh
```

资源默认值位于 `hpc/config/stage_resources.json`。若集群节点内存或 walltime 政策不同，可以调整调度资源，但任何改动都必须保证：`cpus >= max(nprocs, NWorkers)`，且 `nprocs × %maxcore` 留有至少约 10% 的节点内存余量。
