# 精确计算

本目录保存两种目标化合物的 ORCA 6.1.1 高精度 NMR 计算包：

1. F3_005，苯甲酰胺（benzamide）；
2. F3_015，(S)-L-脯氨酰胺（L-prolinamide）。

实际 HPC 工作目录为 `ORCA_61_NMR_ULTRA_2CMPD/`。该英文目录名用于避免部分 Linux 调度器、MPI 与临时目录工具对中文路径处理不一致。

重要边界：本目录中的输入、计算谱和实验外观模拟均不是实验数据。只有完成计算、独立实验比较和验证门控后，某一方法分支才能标为 `VALIDATED`；初始交付状态为 `DRAFT_NOT_VALIDATED`。
