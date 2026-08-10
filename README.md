# DK_Computation — ORCA 6.1.1 NMR workflow

本仓库现只保留一个有效计算包：

- [`精确计算/ORCA_61_NMR_ULTRA_2CMPD/`](精确计算/ORCA_61_NMR_ULTRA_2CMPD/)：苯甲酰胺 `F3_005` 和中性 `(S)`-脯氨酰胺 `F3_015` 的 `¹H/¹³C` ORCA 6.1.1 NMR 工作流。

工作流包含显式溶剂 GOAT 构象搜索、PBE0-D4 优化/频率、QRRHO/Boltzmann 构象平均、同方法 TMS 定标、`revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore` 屏蔽常数生产候选、J 耦合和受限的高层诊断分支。

> **当前状态：`INPUTS_PREPARED / HPC_NOT_RUN / DRAFT_NOT_VALIDATED`**
>
> 仓库中是已准备的输入和调度工具，不是已完成的服务器计算结果。噪声、基线、溶剂峰、峰高和线宽仅属于明确标记的 synthetic appearance 层，不是静态 ORCA 第一性原理观测量。

请从 [计算包 README](精确计算/ORCA_61_NMR_ULTRA_2CMPD/README.md) 开始，服务器操作见 [SERVER_QUICKSTART](精确计算/ORCA_61_NMR_ULTRA_2CMPD/docs/SERVER_QUICKSTART.md)。

ORCA、MPI、xTB、许可证、密码、Token 和 SSH 私钥均不包含在仓库中。
