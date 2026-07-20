# DK_Computation — ORCA Figure 3 HPC workflow

本仓库保存 Figure 3 中 52 种化合物的 ORCA 6.1.1 高质量计算工作流。它面向 Linux 超算环境，可通过 SLURM 或 PBS Pro 逐阶段运行。

> **科学边界：**本仓库生成的 NMR、精确质量、同位素包络和谱图外观均为计算或合成模型，不是实验测试数据。人为加入的溶剂峰、水峰、噪声、基线和线宽不会把计算谱图变成实测谱图。

## 仓库内容

- `ORCA_52_HPC_portable_package/`：完整、可审计的计算文件夹。
- `archives/ORCA_52_HPC_portable_package_20260720.tar.gz`：推荐上传至 Linux 超算的归档。
- `archives/ORCA_52_HPC_portable_package_20260720.zip`：UTF-8 文件名兼容的 ZIP 归档。
- `archives/ORCA_52_HPC_portable_package_20260720_ARCHIVE_SHA256.txt`：两个归档的 SHA-256。

完整包包括：

- 52 种化合物；
- 112 种固定的物种、质子化状态、显式溶剂和连续介质情景；
- 1,120 个 GOAT 构象搜索输入；
- SLURM 与 PBS Pro 作业数组；
- 节点本地 scratch、输入哈希、原子回传、失败保留和断点续算；
- GOAT、DFT 优化/频率、QRRHO/Boltzmann、参考物、shielding、J 耦合、NMR 后处理和公式级 HRMS 分支；
- CSV、TXT、JSON、PNG 和 SVG 输出接口。

ORCA、xTB、MPI、许可证文件和任何私钥均未包含在仓库中。

## 上传后验证

推荐使用 `tar.gz`：

```bash
tar -xzf ORCA_52_HPC_portable_package_20260720.tar.gz
cd ORCA_52_HPC_portable_package

python3 package_metadata/audit_portable_package.py
sha256sum -c package_metadata/PACKAGE_SHA256SUMS.txt

chmod +x hpc/submit/*.sh hpc/bin/*.py hpc/tests/*.sh
PYTHONDONTWRITEBYTECODE=1 bash hpc/tests/run_tests.sh
```

预期结果是 portability audit `PASS`，以及 HPC 离线测试 `8/8 PASS`。

## 超算站点配置

首先编辑：

```text
hpc/config/site_config.json
```

至少填写：

- Linux ORCA 6.1.1 主程序的绝对路径；
- 与 ORCA 发行包兼容的 xTB 路径；
- module 加载命令；
- SLURM account/partition/qos，或 PBS project/queue；
- 可选的节点本地 scratch 路径；
- 如有条件，记录 ORCA 二进制 SHA-256。

Windows 版 ORCA 不能在 Linux 节点运行。ORCA 主程序应通过绝对路径直接调用，不应在外层使用 `mpirun` 或 `mpiexec`。具体要求参见 [ORCA parallel calculations](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html) 和 [ORCA troubleshooting](https://www.faccts.de/docs/orca/6.1/manual/contents/quickstartguide/troubleshooting.html)。

## 先运行两个 GOAT 任务

SLURM：

```bash
export HPC_SCHEDULER=slurm

hpc/submit/submit_goat.sh \
  data/GOAT_explicit_job_manifest.csv \
  --array-range 1-2 \
  --max-concurrent 2
```

PBS Pro 将 `HPC_SCHEDULER` 改为 `pbs`。

上述命令默认只生成调度脚本，不会提交。人工检查账号、队列、CPU、内存、walltime、数组范围、ORCA 路径和 scratch 后，再在同一命令末尾增加：

```bash
--submit
```

建议按以下顺序放大：

1. `1-2`；
2. `1-10`；
3. 删除 `--array-range`，按站点配额执行全量任务。

增大或减小 `--max-concurrent` 只改变同时运行的作业数，不改变单个作业的理论方法和质量。

## 完整阶段顺序

1. 固定结构、立体化学、质子化状态、反离子、pD 和溶剂情景；
2. GOAT 构象搜索；
3. 人工检查连接关系、手性、质子转移、溶剂簇和反离子，并批准构象选择 gate；
4. r2SCAN-3c 优化与频率；
5. 检查虚频、连接关系、QRRHO 和 Boltzmann 权重；
6. 独立完成并批准参考物优化、频率和 shielding；
7. 运行分析物 GIAO shielding；
8. 在代表集验证后运行 J coupling；
9. 分别输出原始 shielding/J、clean calculated spectrum 和 synthetic acquisition appearance；
10. 独立生成分子式、精确质量和同位素包络预测；
11. 执行最终哈希与科学边界审计。

结构和构象依赖的 gate 不能在计算前自动批准。作者对初始结构假设的批准不能替代对计算后连接关系、虚频和构象集合的检查。

## 冻结包 QA

- 1,120 个 GOAT 作业；
- 52 个化合物、112 个情景；
- 6,866 个包内文件；
- 输入、坐标和 bundle SHA-256 全部匹配；
- HPC 离线测试 8/8 通过；
- Python AST 40/40、JSON 1,154/1,154 通过；
- 后处理 dry-run 通过，未启动 ORCA；
- 0 个 Windows 运行时绝对路径；
- 0 个缓存、锁、ORCA 二进制或计算输出混入；
- `tar.gz` 与 ZIP 均与源文件夹 6,866 个文件逐字节一致。

更详细说明见：

- `ORCA_52_HPC_portable_package/README.md`
- `ORCA_52_HPC_portable_package/hpc/README.md`
- `ORCA_52_HPC_portable_package/docs/WORKFLOW_AND_MANUAL_GATES.md`
- `ORCA_52_HPC_portable_package/NOTICE_CALCULATED_ONLY.txt`

