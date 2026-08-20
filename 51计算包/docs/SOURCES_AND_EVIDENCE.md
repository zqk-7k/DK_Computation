# 方法依据与官方文档

本页说明各来源支持什么，以及不支持什么。引用用于设计可审计的工作流，不等同于对 51 个目标化合物的实验验证。

## ORCA 官方文档

1. [ORCA 6.1 NMR tutorial](https://www.faccts.de/docs/orca/6.1/tutorials/spec/NMR.html) 说明 ORCA 的 NMR shielding/J 计算入口与示例。具体关键词仍需在固定的 ORCA 6.1.1 版本上做输入测试。
2. [ORCA 6.1 GOAT manual](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/goat.html) 说明 GOAT 的 worker/global iteration 逻辑。GOAT 的内部停止不替代本包的物种、motif、波次和谱稳定门。
3. [ORCA 6.1 parallel execution](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html) 是使用完整 ORCA driver 路径、匹配 OpenMPI 和本地 scratch 的依据。它明确影响运行方式，但不为具体集群保证最佳 ranks。

## 屏蔽常数和基组

- Schattenberg and Kaupp, *J. Chem. Theory Comput.* (2021), [DOI: 10.1021/acs.jctc.1c00919](https://doi.org/10.1021/acs.jctc.1c00919)。该工作在高水平理论参考上评估核屏蔽方法，支持将 DSD 类方法作为高成本验证候选。它不是复杂溶液实验的普适准确度保证。
- Liang *et al.*, *J. Chem. Theory Comput.*, [DOI: 10.1021/acs.jctc.2c00933](https://doi.org/10.1021/acs.jctc.2c00933)。该工作评估 pcSseg 系列和局部致密/复合基组策略，支持 NMR 专用基组的系统选择与预检。含重元素仍需单独路线。
- DELTA50, *Molecules* (2023), [DOI: 10.3390/molecules28062449](https://doi.org/10.3390/molecules28062449)。该实验位移 benchmark 显示最佳方法会随核种、方法、溶剂和回归策略变化，支持独立校准、盲测和不把单一高成本方法视为普适最优。

## 构象与热力学

- Yesiltepe *et al.*, *J. Cheminformatics* (2018), [DOI: 10.1186/s13321-018-0305-8](https://doi.org/10.1186/s13321-018-0305-8)。支持自动化、构象平均和完整 provenance 的 NMR workflow 思路。具体方法和阈值必须重新验证。
- Grimme, *Chem. Eur. J.* (2012), [DOI: 10.1002/chem.201200497](https://doi.org/10.1002/chem.201200497)。为低频振动的 quasi-rigid-rotor/harmonic-oscillator 处理提供方法依据。本包固定 100 cm⁻¹ 与指数 4 作为工程协议，并要求参数哈希和敏感性报告；它不是唯一可选参数。

## 证据边界

这些来源支持：方法候选、基组选择、构象平均、QRRHO 和运行规范。它们不直接支持：

- 一个显式水等于真实样品含水量；
- benzamide 的表现可直接外推到柔性药物、盐或含碘分子；
- revDSD/NoFrozenCore 一定优于 PBE0 的溶液实验位移；
- 静态 ORCA 可准确产生真实噪声、基线、峰高或交换线宽；
- 目标实验谱可用于逐峰修正后仍称为盲预测。

每个结果必须把文献方法依据与本体系的实际独立验证分开陈述。
