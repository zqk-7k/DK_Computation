# 软件与执行说明

主入口是 `run_supplement.sh`，或软件就绪后执行 `scripts/supplement.py execute`。scripts/workflow.py是复用的FEFF任务管理组件，不是本补算包的主入口。

## Quantum ESPRESSO

先使用平台已有的MPI版本。没有可用版本时，install_qe.sh从QEF官方仓库下载7.4.1，检查固定commit：

`500de340b820e1cb8c05f2d8bb8fced102f377c1`

官方入口：<https://github.com/QEF/q-e/releases/tag/qe-7.4.1>；构建说明：<https://www.quantum-espresso.org/Doc/user_guide/node11.html>。

需要匹配的C/Fortran/MPI编译器、Git、CMake≥3.20和make。默认只编译pw目标，启用MPI、关闭OpenMP，并使用内部FFT/LAPACK作为可移植基准。平台已有优化BLAS/FFT时，可由Codex按平台文档配置，记录CMakeCache和修改，不能混用不兼容MPI。编译并行数由XAS_BUILD_JOBS控制，默认2；遵守登录节点编译政策。

已有pw.x的MPI帮助横幅会在自动启动前检查。若无法确认MPI版本，流程停止并要求检查模块，不会将串行程序重复启动16份到同一scratch目录。用户显式配置dft.command时应确保命令与所申请MPI资源匹配。

当前DFT方法不是从旧CIF还原出的原始VASP方法。新模型默认PBE自旋极化/U0；U、磁态、溶剂、表面质子化和电位对应关系需要后续方法敏感性检查。QE卡片说明：<https://www.quantum-espresso.org/Doc/INPUT_PW.html>。

## 赝势

使用官方Materials Cloud存档的SSSP1.3.0 PBE efficiency。版本固定，不自动换成最新库：

<https://archive.materialscloud.org/records/rcyfm-68h65>

元数据MD5：`3153c4b20fc90a44fba0236627525644`；赝势压缩包MD5：`a58f1b3373f330179fd0832c48bb9a52`。安装程序只提取H/C/O/Co/Fe/W六种UPF，并逐文件核对官方摘要，另记录SHA256。不会解压外部路径或运行压缩包中的代码。

所选Fe赝势的建议波函数/电荷密度截断为90/1080Ry；因此统一采用不低于此值的基准。建议值不代替具体体系的收敛验证。若替换赝势库，所有比较模型需要使用一致且经验证的设置，保留原结果。

## FEFF10与Python

FEFF10沿用上一包的官方固定commit和SHA256，安装脚本中完整记录。优先复用已有可工作的FEFF10模块集合；软件命令必须依次完成相应模式所需模块，不能仅指向单个模块。官方仓库：<https://github.com/times-software/feff10>；手册：<https://times-software.github.io/feff10/assets/doc/feff10_users_guide.pdf>。

Python采用已运行过的无界面数值依赖路线：固定NumPy/SciPy/Matplotlib、pymatgen和xraylarch。xraylarch以--no-deps安装，明确省略GUI/XRD依赖；pip check的缺项写入日志，Larch数值导入及真实xftf/cauchy_wavelet自检必须通过。完整环境版本冻结到software/python_environment.lock.txt。WT实现来自Larch：<https://xraypy.github.io/xraylarch/xafs_wavelets.html>。

## 平台适配和续跑

config.json中的队列和账户默认空值，表示采用平台默认。不要猜造队列；先查询实际系统。上一轮使用过Slurm，当前仍由doctor重新确认。module_commands用于计算节点加载真实模块，安装阶段则先在当前会话加载同一套模块。

DFT默认48小时队列时限，程序超时约46小时；FEFF默认72小时队列时限，程序超时70小时。超时只是资源限制，不能充当预计耗时。新XANES网格较宽，应先检查代表任务的实际运行情况。

计算中断后可重跑主命令；程序依据输入、输出和结构摘要跳过已完成任务，并先查询已有数组。任务锁存在时先检查进程/调度器，不可盲删锁后重复提交。失败的优化不能进入“已优化”结构集；正常退出但未满足SCF/离子/力收敛也会被拒绝。

若改变DFT参数或结构，应复制成新计算目录，以便保留旧配置、日志和结果。主流程对这些更改采用拒绝覆盖策略。调度器暂时不可查询时，恢复查询后再续跑，不重复提交未知状态的任务。

## 成功的定义

程序验证与物理计算分开记录。本地自检只验证代码、数据格式和图件；只有实际完成QE优化、FEFF计算、宽能区检查、位点平均和完整导出后，才生成标为完成的结果包。结果报告还会单独说明是否保留了目标O配位，不能把“程序完成”解释成“目标机理已证明”。
