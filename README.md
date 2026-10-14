# DK_Computation — ORCA NMR workflows

当前 64 个化合物项目：

- [`工作2-计算资源包-5/`](工作2-计算资源包-5/)：64 个化合物、66 个结构分支，沿用 SCNET ORCA 6.1.0 wrapper。入口为 [先读说明](工作2-计算资源包-5/README_先读.txt)，本次调度、恢复和默认 30 天墙时修改见 [生产推进与长墙时修订](工作2-计算资源包-5/docs/生产推进与长墙时修订_20261014.txt)。
- 该项目当前完成的是代码与本地测试修订；SCNET 真实量化小试和 64 化合物生产计算尚未由本轮启动。提交前必须完成真实小试，默认墙时也必须通过实时调度校验。

历史两化合物项目（独立保留）：

- [`精确计算/ORCA_61_NMR_ULTRA_2CMPD/`](精确计算/ORCA_61_NMR_ULTRA_2CMPD/)：苯甲酰胺 `F3_005` 和中性 `(S)`-脯氨酰胺 `F3_015` 的 `¹H/¹³C` ORCA 6.1.1 NMR 工作流。

历史两化合物工作流包含显式溶剂 GOAT 构象搜索、PBE0-D4 优化/频率、QRRHO/Boltzmann 构象平均、同方法 TMS 定标、`revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore` 屏蔽常数生产候选、J 耦合和受限的高层诊断分支。

> **当前状态：`INPUTS_PREPARED / HPC_NOT_RUN / DRAFT_NOT_VALIDATED`**
>
> 仓库中是已准备的输入和调度工具，不是已完成的服务器计算结果。噪声、基线、溶剂峰、峰高和线宽仅属于明确标记的 synthetic appearance 层，不是静态 ORCA 第一性原理观测量。

历史项目入口为 [计算包 README](精确计算/ORCA_61_NMR_ULTRA_2CMPD/README.md)，服务器操作见 [SERVER_QUICKSTART](精确计算/ORCA_61_NMR_ULTRA_2CMPD/docs/SERVER_QUICKSTART.md)。两个项目的软件版本与配置不可混用。

ORCA、MPI、xTB、许可证、密码、Token 和 SSH 私钥均不包含在仓库中。
