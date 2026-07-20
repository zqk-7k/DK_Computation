# 52 个化合物：24 小时隐性溶剂 NMR/HRMS 计算包

## 这套程序能完成什么

这是一套可直接上传到 GitHub、再克隆到 Ubuntu 超算节点运行的完整程序。它面向目前确认的 52 个 Figure 3 结构，目标是在一次独占节点 allocation 开始后的 24 小时内完成：

- 52 个化合物的公式级高分辨质谱精确质量和天然同位素包络；
- 所有具有可观察氢的化合物的计算型 ^1H NMR；
- 所有含碳化合物的计算型 ^13C{^1H} NMR；
- 每个谱图的完整点阵 CSV、TXT、PNG、SVG 和 Excel；
- 原始 ORCA 输出、计算方法层级、结构哈希、参考物结果和最终审计表。

这里的“完整谱图”表示：ORCA 计算化学位移，几何/拓扑经验模型生成限时条件下的 ^1H 裂分，随后生成干净计算谱和带噪声、基线、溶剂峰及水峰的合成外观谱。经验 J 并不是 ORCA 全矩阵从头算 J。所有文件都会标为 `CALCULATED/SYNTHETIC - NOT EXPERIMENTAL`。

HRMS 不调用 ORCA。它按离子分子式、精确同位素质量和天然丰度卷积生成同位素包络，不模拟碎裂、离子丰度、污染峰或仪器实测噪声。

## 为什么能把时间控制在 24 小时内

程序采用“先全覆盖、再升级”的顺序：

1. 用 xTB/ALPB 隐性溶剂优化每个溶质结构。
2. 先对全部可观察物种运行 `PBE0/pcSseg-1/CPCM` shielding。
3. 若主方法失败，依次使用明确标记的 `PBE/pcSseg-1` 和 `PBE/def2-SVP` 保底方法。
4. 全覆盖完成后，按大分子和柔性分子优先，运行限时 GOAT、`r2SCAN-3c/CPCM` 优化和 `PBE0/pcSseg-2/CPCM` shielding。
5. 4 小时为几何门限，6 小时为参比门限，14 小时为主基线门限，16 小时为紧急保底门限；19 小时后停止启动新的精修，21.5 小时终止仍未完成的可选精修；23.25 小时前完成谱图和审计，23.5 小时为程序硬截止，Slurm 在 23 小时 55 分钟终止作业。

最终结果会逐个标明 `TIER1_DFT_SHIFT`、`TIER2_DFT_SHIFT` 或 fallback。未完成的精修不会覆盖已完成的基线，也不会被误标成高精度结果。

“24 小时”从节点真正开始运行本程序时计时，不包括 Slurm 排队时间。若必须从提交命令开始计时，需要预留节点。

## 固定的资源配置

当前配置针对：

- 2 × Intel Xeon Platinum 8470Q；
- 104 个物理核、208 个逻辑线程；
- 754 GiB 内存；
- 本地 NVMe；
- Ubuntu 22.04.5。

程序最多使用 92 个物理核，保留 12 个核给系统和 I/O，低于整机物理核的 90%。它禁用超线程计算，并将 OMP、MKL、OpenBLAS 和 BLIS 固定为单线程。RTX 5090 不参与当前 ORCA/xTB/DFT 路线。

## 科学模型边界

- 全部电子结构计算只使用隐性溶剂：DMSO-d6 对应 `DMSO`，D2O 对应 `WATER`。
- F3_016 和 F3_017 使用去除旁观 Cl− 的 +1 阳离子连续介质模型。
- F3_021 使用去除两个旁观 Cl− 的 +2 二阳离子连续介质模型。
- F3_019、F3_020 和 F3_023 使用已确认的两性离子。
- F3_035、F3_043 和 F3_044 不含碳，因此 ^13C 输出为 `NOT_APPLICABLE_NO_CARBON`，不会伪造碳峰。
- F3_035、F3_043 和 F3_044 在 D2O 中没有持续可观察的分析物氢；^1H 合成外观只保留 HOD/背景，不生成虚假分析物峰。
- F3_006 含 I。pcSseg 只到 Kr，因此 I 使用 `def2-TZVP/def2-ECP`，而目标 H/C 仍使用 pcSseg。
- D2O 中的 NH/OH 交换氢不会进入分析物 ^1H 峰表。
- DMSO 中 NH/OH 的位置对特异氢键和浓度敏感。纯隐性溶剂只能给出低置信度预测，程序会保留宽峰并明确来源。
- 24 小时版本每个化合物只保留一个基线构象或一个 GOAT 最低能构象；它不是构象系综的 Boltzmann 平均。柔性分子和酰胺旋转异构体因此具有更高不确定性。
- ^1H 裂分使用限时的一级经验 J 模型，不是完整自旋 Hamiltonian。强耦合、二级谱效应和实验浓度引起的交换行为不能由本流程精确再现。

## 服务器上一次性运行

### 1. 克隆仓库并进入计算包

```bash
git clone git@github.com:zqk-7k/DK_Computation.git
cd DK_Computation/ORCA_52_IMPLICIT_24H
```

### 2. 建立 Python 环境

```bash
bash hpc/bootstrap_environment.sh
```

这一步需要能访问 Python 包源。若计算节点不能联网，请在登录节点执行，或让管理员提供 RDKit、NumPy、SciPy、Matplotlib、openpyxl 和 psutil。

### 3. 确认 ORCA

需要 ORCA 6.1.1 Linux x86-64 shared OpenMPI 4.1.6 版本及与它匹配的 OpenMPI。请使用 ORCA 主程序的绝对路径，例如：

```bash
/opt/orca_6_1_1/orca
```

不要在外层运行 `mpirun orca`。ORCA 会根据输入中的 `%pal nprocs ... end` 自行启动并行模块。

### 4. 提交 Slurm 作业

```bash
bash hpc/submit_24h.sh /opt/orca_6_1_1/orca
```

若集群要求分区、账户或 QOS，可直接追加 sbatch 参数：

```bash
bash hpc/submit_24h.sh /opt/orca_6_1_1/orca \
  --partition=cpu \
  --account=YOUR_ACCOUNT \
  --qos=normal
```

提交脚本申请单节点、92 个物理核、680 GiB 内存、23 小时 55 分钟，并在程序内部执行严格时间门控。生产作业要求 `SLURM_TMPDIR`，或者预先导出位于节点本地 NVMe 的 `NVME_SCRATCH_ROOT`；预检要求至少 300 GiB 可用空间。

### 5. 查看状态

```bash
squeue -u "$USER"
tail -f logs/slurm-JOBID.out
bash hpc/status.sh results/production_JOBID
```

主要进度文件是：

```text
results/production_JOBID/pipeline_status.json
results/production_JOBID/baseline_coverage.json
results/production_JOBID/pipeline_final.json
results/production_JOBID/FINAL_AUDIT.json
results/production_JOBID/FINAL_COMPOUND_AUDIT.csv
results/production_JOBID/FINAL_FILE_MANIFEST_SHA256.csv
```

### 6. 打包结果

作业完成后运行：

```bash
bash hpc/make_result_archive.sh results/production_JOBID
```

它会生成 `.tar.gz` 和对应的 SHA-256。

## 没有 Slurm 时

确保获得独占节点后运行：

```bash
export NVME_SCRATCH_ROOT=/local_nvme/$USER
bash hpc/run_direct_24h.sh /opt/orca_6_1_1/orca
```

## 小规模检查

不调用 ORCA，只检查结构和输入模板：

```bash
python scripts/run_24h_pipeline.py \
  --dry-run \
  --only F3_001,F3_006,F3_021,F3_022 \
  --run-dir tests/dry_run
```

在小节点做真实 pilot：

```bash
python scripts/run_24h_pipeline.py \
  --orca /opt/orca_6_1_1/orca \
  --scratch /local_nvme/orca_pilot \
  --run-dir results/pilot \
  --only F3_001,F3_006,F3_021,F3_022 \
  --core-cap 16 \
  --deadline-hours 4
```

## 输出结构

```text
results/production_JOBID/
├── HRMS_CALCULATED_NOT_EXPERIMENTAL/
├── NMR_CALCULATED_NOT_EXPERIMENTAL/
│   ├── F3_001/
│   │   ├── 1H_predicted_peaks.csv
│   │   ├── 1H_spectrum_full_points.csv
│   │   ├── 1H_spectrum_full_points.txt
│   │   ├── 1H_clean.png / .svg
│   │   ├── 1H_synthetic_appearance.png / .svg
│   │   ├── 13C_predicted_peaks.csv
│   │   ├── 13C_spectrum_full_points.csv
│   │   ├── 13C_spectrum_full_points.txt
│   │   ├── 13C_clean.png / .svg
│   │   ├── 13C_synthetic_appearance.png / .svg
│   │   ├── NMR_prediction_data.xlsx
│   │   └── metadata.json
│   └── ... F3_052/
├── raw_orca/
├── parsed_shieldings/
├── geometries/
├── FINAL_COMPOUND_AUDIT.csv
├── FINAL_AUDIT.json
└── FINAL_FILE_MANIFEST_SHA256.csv
```

## 官方方法依据

- ORCA NMR 教程：<https://www.faccts.de/docs/orca/6.1/tutorials/spec/NMR.html>
- ORCA GOAT 手册：<https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/goat.html>
- ORCA 隐性溶剂手册：<https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html>
- ORCA 并行运行手册：<https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html>

## 不能承诺的部分

程序可以严格控制资源和截止行为，但任何电子结构计算都可能遇到异常 SCF、几何变化、文件系统故障或节点中断。因此不能在运行目标服务器 pilot 前对“52/52 都成功”作数学保证，也不能保证 52 个分子都达到 Tier 2。程序的硬目标是先得到完整 Tier 1/fallback 覆盖，再尽可能多地升级。最终 `FINAL_AUDIT.json` 为 `PASS`、`pipeline_final.json` 为 `COMPLETE`，且最终 SHA-256 清单存在，才表示该次结果完整。

计算谱可用于实验规划、方法比较和大模型评估。正式投稿时必须用真实仪器数据替换，不能把这些预测谱称为实测表征。
