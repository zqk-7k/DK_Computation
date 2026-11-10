# W–O 配位与三状态 XAS 补算流程包

本包可交给 Linux 超算端 Codex 执行，包含候选结构构建、Quantum ESPRESSO 优化、FEFF10 计算、Larch 后处理、TXT 和图片导出。原始16个CIF作为不可改写的参考保存；另行构建12个待优化候选模型。

**这里交付的是流程和未优化候选结构，不是已完成的补算结果。** 新几何必须经过实际DFT优化；光谱必须由实际FEFF输出产生。所有输出按计算结构命名，不预设“Co不变化、W变化”的结果。

## 对三个要求的落实

1. **W 第一配位壳为 O：** 原来的顶吸附W结构含有较近的Co/Fe邻居，不能直接改谱。本包从原始CoFe LDH另建W替换Fe或Co的两个候选位点。每个初始结构的W第一重原子壳为六个O，优化中全部原子自由移动，晶胞固定。优化后重新检查真实近邻和键长。
2. **W L₃ 能区10100–10300 eV和W–O FT峰：** 使用FEFF自身的宽能量网格计算，并强制检查原生输出覆盖范围。缺失能区不会用外推、补零或移动单条谱线填补。输出真实W–O几何距离和未作相移校正的FT峰位置；约1.5 Å仅作为参考图的观察区域，不作为优化约束。
3. **W、Co、Fe三行对比图：** 每行叠加羟基化、Co邻位去质子化、EG结合三种结构的计算谱。另保留EG结合Co、EG加入去质子化背景等对照。默认导出8组对比，避免只选支持某一活性位点结论的结构。

尚无CHE/恒电位或实验结构证据把这些候选映射到OCP和1.6 V，所以图例使用结构状态名称。申请的电化学条件保留在config.json中作为待验证对应关系。不能仅凭某元素谱线变化，证明该元素是唯一活性位点。

## 最简单的使用方法

把解压后的整个目录上传到超算，建议改用短的英文目录名，例如 `~/projects/xas_o_coord`。把 `CODEX_TASK.txt` 发给超算端Codex，让它先读取 `AGENTS.md`。

```bash
bash run_supplement.sh auto
```

Codex应先确认平台模块、账户和队列，写入config.json；已有合适软件时优先使用。程序支持Slurm和OpenPBS/PBS Pro，Torque数组语法也可配置。实际DFT/FEFF只在计算节点运行；找不到调度系统时不会自动在登录节点启动计算。

首次安装使用项目目录中的独立Python环境，固定FEFF10版本、Quantum ESPRESSO 7.4.1和SSSP 1.3.0 PBE efficiency赝势。安装脚本不需要sudo。QE优先使用平台MPI版本；默认多核运行前检查MPI标识，避免把串行程序重复启动到同一工作目录。

软件获取和安装说明见 `docs/SOFTWARE_AND_EXECUTION.md`。离线平台可从允许外连的节点下载，再把文件放入software目录；不得绕过平台访问策略。

## 顺序与结果

| 步骤 | 执行内容 | 主要文件 |
|---|---|---|
| 1 | 构建、核对12个配对候选 | candidates/*/structure.cif、metadata.json |
| 2 | 固定晶胞、原子位置DFT松弛 | dft/*/qe.in、qe.out、relaxed.cif |
| 3 | SCF/离子收敛、最终力和W配位检查 | result.json、validation/relaxed_coordination.json |
| 4 | 两个代表模型的FEFF半径检查 | screen_small、screen_large，共24个任务 |
| 5 | 对每个模型所有W/Co/Fe吸收位点计算 | supplement，共384个XANES/EXAFS任务 |
| 6 | 位点平均、FT和WT、三行图 | results/supplement/COMPLETE_*/ |
| 7 | TXT、图片和完整结果包 | FINAL_EXPORT_*、deliverables/*.zip |

`FINAL_EXPORT_*/txt/`包含全部位点和平均谱的TXT表格、三状态同轴数据、FT实部/虚部/幅值、WT长表、差谱指标及W第一壳峰位。单位和来源写在表头，Origin可直接导入。

`FINAL_EXPORT_*/figures/`包含PNG、PDF、SVG、TIFF。图形为白底、统一字体，优先Arial；平台没有Arial时如实记录替代字体。W/Co/Fe按三行排列，颜色在各图保持一致。

优化后如果有模型不再保持纯O第一壳，默认仍保留该结构并计算其谱，配位审计和完成报告会标明假设未保持。不能把这些模型称为已确认的W–O结构。若希望在此处停止，可将 `coordination.failure_policy` 改为 `stop`。两种模式都不修改原子距离来迎合目标峰位。

## 计算方法与资源

对于非对称slab，默认自动把零外场偶极修正窗口放入真空区，并在优化后复核原子未进入该区域；这项修正不代表施加1.6 V电极电位。仍需检验真空厚度和静电势平台。

当前新DFT基准为自旋极化PBE、U=0、固定晶胞、全部原子松弛、2×2×1 k网格、90/1080 Ry截断。这是明确的新计算设置，不是从原CIF恢复出的VASP参数。90/1080 Ry满足所选SSSP中Fe的建议值，但仍需做本体系的截断、k网格、磁态及U敏感性检验。当前流程没有溶剂、电极电位或工况下结构占比模型，不能自动给出1.6 V的稳态结构。详见 `docs/STRUCTURE_AND_METHOD_LIMITS.md`。

默认资源：DFT每任务16个MPI进程、64 GB、最多2个并行任务；FEFF每任务1核、4 GB、最多16个并行任务。可在确认配额后调整。更多FEFF任务并行比给串行FEFF单任务增加核数更有效。

上一轮真实结果中，XANES任务中位耗时约12小时、EXAFS约10小时；这取代了此前未实测的30–60分钟假设。此次还有DFT优化、更宽能量网格和较大FMS半径，不能承诺几天完成。必须先看代表任务的实际耗时、内存与收敛，再调整72小时的FEFF任务时限和并发；总体可能需要数天至数周，排队另计。若超时，保留结果并在新设置下明确重试，不能将部分输出标为完成。

FT/WT是对同一chi(k)的数学处理，不是额外独立测量。FEFF从几何计算自身电子势；本流程未把QE电荷密度转入FEFF。解释价态和活性位点时必须同时考虑这个方法边界。

## 分步命令

```bash
# 准备候选，便于检查；不进行物理计算
.venv/bin/python scripts/supplement.py prepare --config config.json
# 检查环境与原始文件摘要
.venv/bin/python scripts/supplement.py doctor --config config.json
# 软件已安装后，运行和续跑完整流程
.venv/bin/python scripts/supplement.py execute --config config.json --scheduler slurm
# 已有完整后处理结果时单独重出TXT和图
.venv/bin/python scripts/export_supplement.py --config config.json --result-dir results/supplement/COMPLETE_实际时间
```

长任务监控可放在平台允许的tmux会话中。不同U、磁态、赝势、数值设置必须复制为新计算目录，保留旧输入和结果；不可直接覆盖正在运行的任务。流程会检查输入和输出摘要，避免错误复用结果。

## 本地验证范围

本机已检查候选构建、原子身份与化学计量、QE输入/输出解析、FEFF宽能区输入、结果完整性检查、TXT导出和图件。测试用谱只存在临时目录，不是交付的计算结果。

本机没有运行新DFT优化或新FEFF物理计算，也没有验证目标超算的MPI/队列环境。实际运行从上传后的软件与环境检查开始。可追溯验证记录见validation目录。
