# CIF → FEFF → XAS 超算流程包

本包用于由已优化的原子结构计算理论 XAS。输入的16个CIF已原样保存并记录校验值。用户已确认结构完成Materials Studio/VASP优化，因此默认不重复DFT。

**本包交付的是可执行流程与输入，不是已经完成的FEFF计算结果。** 原DFT的U值、磁性、电子数和电荷密度未提供。FEFF将根据几何结构及本包明确的设置计算自身局域电子结构；不能声称继承了原VASP电子态或外加电位。

## 先上传什么

将整个解压后的目录上传到Linux超算可读写的项目目录，建议使用短的英文路径，例如 `~/projects/w_cofe_xas`。保留 `structures/`、`scripts/`、`config.json` 等目录的相对关系。

把 `CODEX_TASK.txt` 的内容发给超算上的Codex。`AGENTS.md`已经包含执行规则，Codex应先读取。

## 推荐执行顺序

1. 检查平台：识别Slurm/PBS、默认队列、账户、内存、时间限制、编译器模块和网络条件。
2. 用 `config.json` 记录实际队列与资源；空队列/账户表示采用平台默认值，不是虚构队列。
3. 安装项目独立Python环境和固定版本FEFF10，不需要sudo，不改全局环境。
4. 准备代表位点的6个任务：W L3、Co K、Fe K各一位点，各计算XANES和EXAFS。
5. 比较增大SCF/FMS/原子团半径前后的实际输出，按配置阈值判断数值变化。通过后采用较大的半径进入生产任务。
6. 默认计算核心7个结构，每个结构内所有Co/Fe/W吸收原子逐个计算；共232个小任务。串行FEFF任务由调度数组并行，默认最多同时4个。
7. 验证真实输出、完成位点平均、生成XANES、加权EXAFS、复数FT和WT。
8. 检查图件和输出清单；报告收敛范围、方法假设及仍需实验检验的问题。

自动主流程命令：

```bash
bash run_all.sh core auto
```

它会安装依赖、检查环境、进行后处理自检，然后提交并等待实际计算完成。未知队列/账户导致提交失败时，Codex应读取错误、查询平台并修正配置，不能猜造设置。网络禁止外连时使用平台允许的下载节点/预置模块；不能声称下载已成功。

长时间运行建议让协调程序保持在 `tmux` 会话内；监控程序本身很轻，实际FEFF计算必须进入队列。断线后重新执行同一命令，会依据真实输出和输入哈希跳过已完成任务，并继续监控已有数组。

## 为什么先做7个，而不是16个

核心组为：

| 结构文件 | 作用 |
|---|---|
| CoFe-LDH.cif | 无W基底 |
| Top-O1.cif | 一种W锚定环境，近邻包含Co/Fe |
| Top-O2.cif | 另一种W锚定环境，用于配位环境对照 |
| W1-path/CH2OH-CH2OH.cif | W体系的EG吸附结构 |
| W1-path/CH2OH-CHO.cif | W体系的关键羟基乙醛中间体结构 |
| CoFe-LDH-path/CH2OH-CH2OH.cif | 无W体系的EG吸附对照 |
| CoFe-LDH-path/CH2OH-CHO.cif | 无W体系的关键中间体对照 |

剩余9个保留为后续反应路径扩展，只有研究后续中间体或核心模型不能解释谱图时才需要全部计算。`2HCOOH+sub.cif`更接近产物离开表面的构型，初筛优先级较低。

全部16个结构可用：

```bash
bash run_all.sh full auto
```

full为526个任务。不要把core和full同时无差别提交；先看核心结果是否需要扩展。

## 常用分步命令

```bash
bash scripts/bootstrap_python.sh
bash scripts/install_feff10.sh
.venv/bin/python scripts/workflow.py doctor --config config.json
.venv/bin/python scripts/workflow.py prepare --profile pilot --config config.json
.venv/bin/python scripts/workflow.py run --profile pilot --scheduler slurm --dry-run
.venv/bin/python scripts/workflow.py execute --profile core --scheduler slurm
```

PBS使用 `--scheduler pbs`。默认按OpenPBS/PBS Pro生成 `-J` 数组；Torque应在配置中设置 `pbs_flavor="torque"`。若集群数组语法有定制，由Codex根据平台文档调整并先dry-run。

只有在已经申请到计算节点的交互会话内，才可显式本地运行：

```bash
.venv/bin/python scripts/workflow.py execute --profile core --scheduler local --allow-local
```

不能在登录节点执行这条命令。主流程不会在找不到调度器时偷偷改成本地大批计算。

## 配置中需要平台现场确认的项目

- `scheduler.type`：auto/slurm/pbs。
- `partition`、`queue`、`account`：实际平台要求的值，未要求时留空。
- `memory_gb`、`walltime`、`max_concurrent`：默认8GB、4小时/任务、4并发，可根据小任务实测调整。
- `module_commands`：计算节点运行前需要加载的模块命令；只填写平台真实支持的命令。
- 已有FEFF时，可将 `software.feff_command` 设为命令数组，例如 `['/path/to/run_feff10.sh']`（JSON中使用双引号）。该命令必须完成对应任务的全部FEFF模块，不能只填写其中一个模块。

队列信息未知不会影响准备流程包。Codex应首先运行doctor并阅读平台使用说明，查询到缺失参数后继续，而不是反复向用户请求普通安装和检查许可。

## 已明确采用的计算约定

- 几何结构按用户提供CIF保留；不删除H，不压掉真空，不将slab强行转成块体。
- Core/full按同一模型同一元素所有显式位点等权；pilot只用代表位点，不能代表整体Co/Fe信号。
- 初始设置：原子团8Å、SCF5Å、FMS6Å、RPATH5Å；XANES与EXAFS独立输入目录；EXAFS算至k=14Å^-1。
- 生产半径由pilot增大测试选择；其他物理设置（芯空穴、无序、磁性等）仍需科学审阅，数值收敛不等于模型正确。
- 当前是结构基准计算；默认芯空穴FSR，S0²=1作为理论振幅设置，不是从实验拟合出的数值；未引入人为噪声、随机峰或人为电位响应。
- 当前采用未指定自旋极化的基准输入，没有继承VASP中的Co/Fe磁有序。若谱线解释依赖磁性或特定电荷态，应另建有明确设置的对照计算。
- XANES在FEFF绝对光子能量轴上平均，不逐条对白线峰作平移对齐。
- EXAFS在各位点原生阈值下的共同k网格上先平均带符号chi，再FT/WT；这是结构平均近似，未额外构造多位点绝对能量化学位移混合。
- FT保留幅值/实部/虚部；WT来自相同计算chi，不是独立实验证据。
- 理论谱按结构/吸附态命名。现有CIF未建立OCP/1.4/1.5/1.6V对应关系，不使用这些电位标签。
- 未另行加入由实验、声子或分子动力学确定的热无序。与实验振幅作定量比较前，应再设置有依据的温度/无序模型。
- 没有实测谱时，不生成虚构的实验拟合误差或配位数置信区间。

## 会得到什么

- `jobs/`：每个位点每个模式的feff.inp、job.json、原生输出、日志、运行结果。
- `manifests/`：结构、位点、权重和输入的追溯清单。
- `validation/`：环境、收敛筛查、任务完整性和最终运行状态。
- `results/`：各结构三条边（存在该元素时）的XANES、EXAFS、k加权、FT、WT数据及图件。
- CSV/NPZ保留完整数值；PDF/SVG/PNG用于查看及排版。标题标注Calculated FEFF，不把理论输出写成实测。
- `deliverables/CALCULATED_XAS_*.zip`：最终自动打包，含原生文本输出、结构、输入、处理结果和校验清单；可重建的大型临时二进制文件不重复打入ZIP。

原目标五页图集的数据类型均被覆盖；可按真实结果重新排版。不同模型的线形与变化幅度由计算得到，不以旧示意图的趋势为目标函数。

## 软件和版本

FEFF10使用官方固定commit和压缩包SHA256，见 `docs/SOFTWARE_SETUP.txt`。Python核心版本固定于requirements.txt；完整安装环境会保存到software/python_environment.lock.txt。

Athena不是纯理论数据后处理的必需项。FDMNES是后续全势XANES交叉验证的可选项，未混入主FEFF计算，也未假称已经运行。其官方入口及部署说明见docs目录。

## 耗时如何估计

目前不能仅凭CIF准确给出总耗时。先运行6个pilot任务，读取每个`run_result.json`的`duration_seconds`和调度器内存统计，再调整资源。核心组232个任务、全量组526个任务；总CPU时间约为各类任务实测平均时间乘以对应任务数。墙钟时间还取决于并发数、队列等待和增大半径后的开销，不能简单把pilot时间当作生产任务时间。软件安装、半径筛查、重试和后处理需另计。

更改FEFF版本、芯空穴或其他物理设置时，应建立新profile/计算目录，保留原输出与来源，不复用旧的完成状态。

## 本地验证范围

已在独立Python环境检查结构解析、输入生成、原生表格解析、Larch后处理和图件导出。测试用数据只用于程序验证，不作为计算结果交付。

当前机器未提供可运行的Linux FEFF编译环境，因此FEFF10的目标平台编译、实际光谱和超算提交必须在上传后完成。`validation/LOCAL_VALIDATION.json`记录实际完成与未完成的检查。主流程会校验输出而不是仅依赖退出码。
