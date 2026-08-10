# 方法选择、精度预期与限制 / Method rationale and limitations

## 1. “最高配置”不是单一关键词

ORCA 可以运行很昂贵的方法，但 NMR 溶液预测没有一个对所有分子、所有核种和所有条件都最准确的通用“最高档”。总误差同时来自：

- 电子结构方法和基组；
- 构象遗漏与自由能排序；
- 显式/连续溶剂模型；
- 温度、浓度和聚集；
- 质子化、盐型、互变异构和交换；
- 参考物及参考计算一致性；
- 振动、相对论和局域相关近似；
- 实验采集和处理。

因此，本包采用“生产主线 + 全构象交叉检查 + 小子集高阶诊断 + 留出验证”，而不是把最贵的单点计算直接称为最准确结果。

## 2. 构象层

### GOAT-ENTROPY/XTB/ALPB(DMSO)

GOAT 用于广泛搜索低能构象和氢键微结构。每个固定计量场景使用 20 个独立、确定性起点，以降低单一起点遗漏风险。XTB 能承受较大的构象任务量，但其能量排序不是最终自由能，因此所有保留构象必须进入更高层级优化/频率。

### PBE0-D4/def2-TZVPP 最终优化与频率

PBE0-D4 提供杂化泛函、色散和三重 ζ 几何基础。`VeryTightOpt`、`VeryTightSCF` 和 `DefGrid3` 用于减少几何和数值噪声。频率用于确认驻点并提供热化学/QRRHO 输入。它仍不能证明显式溶剂计量就是实验中的真实占比。

## 3. shielding 主线

### revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore

这是本包的生产候选，而不是预先宣布的获胜者。选择理由是：

- 双杂化 DFT 同时包含较高比例交换和 MP2 型相关；
- `pcSseg-3` 针对磁屏蔽性质设计；
- `NoFrozenCore` 避免核心相关处理在不同分支间不一致；
- GIAO 提供规范原点不变的磁响应；
- CPCM(DMSO) 与显式 DMSO 簇共同描述体相和局部氢键。

代价是计算昂贵、显式簇更难收敛，并可能对 MP2 密度、基组和微溶剂结构敏感。它必须与同方法 TMS 参考配对，并由实验校准/留出结果决定是否晋级。

### PBE0/pcSseg-3 全构象交叉检查

PBE0 分支计算成本较低、数值通常更稳。它对所有批准构象运行，用于：

- 发现 revDSD 异常值或收敛问题；
- 估计主要方法依赖性；
- 检查峰顺序是否稳健；
- 在高阶分支失败时提供透明的替代结果。

报告必须同时显示 revDSD–PBE0 的位移差异，不能只挑与实验更接近的一支。

## 4. MP2 与 CCSD(T) 受限基准

### RI-MP2/pcSseg-3

RI-MP2 是小规模 canonical MP2 响应基准，限制在每个固定计量的前 3 个构象。它用于判断双杂化中 MP2 型相关的方向，不承担全构象生产工作。

### TightPNO DLPNO-MP2/pcSseg-3

DLPNO-MP2 以局域相关降低成本，并使用更严格 PNO 设置。它限制在每个固定计量的顶级构象。局域近似本身会引入误差，不能仅因“MP2”名称而假定优于 DFT。

### AUTOCI-CCSD(T)/pcSseg-4

该分支是最高电子结构诊断之一，但仅用于资源门控后的溶质单体/截断显式模型。当前包不声称已经验证 CCSD(T) 的 CPCM NMR 响应，因此把它作为气相或显式截断的 **delta diagnostic**，不直接替代 DMSO 溶液生产值。

CCSD(T) 对 NMR 不会自动消除构象、溶剂、参考和交换误差。对小分子单体极高精度的电子相关，不等于对真实溶液整谱极高精度。

## 5. VPT2 振动修正

VPT2 可估计非谐振动对平均性质的修正，但成本高且容易受低频、共振和显式溶剂软模影响。因此默认仅考虑：

- 每个物种最重要的 1–3 个构象；
- 溶质单体或最多 1 个显式 DMSO 的微结构；
- 完整的数值稳定性与共振审查。

VPT2 不描述 NH 在 DMSO/水中的化学交换动力学，也不保证改善 NH 化学位移或峰宽。

## 6. 显式溶剂与跨计量禁混

苯甲酰胺和脯氨酰胺的羰基、酰胺 NH₂ 以及脯氨酰胺环 NH 对氢键高度敏感。连续介质不能表示特定氢键，因此设置 0–3 DMSO 与痕量水敏感性分支。

不同显式溶剂计量含有不同粒子数，不能用普通封闭体系 Boltzmann 因子直接相加。正确做法是：

- 在每个固定计量内部做构象平均；
- 将不同计量作为并列模型报告；
- 只有在另行建立溶剂化学势、标准态和水活度模型后，才可能推导跨计量总体权重。

`2 DMSO + 1 H₂O` 只回答“少量水可能怎样改变结果”，不预测样品的真实水含量或 H₂O/HDO/D₂O 同位素分布。

## 7. 参考定标

屏蔽常数不是直接可比较的实验化学位移。同一分支必须使用同方法、同基组、同溶剂、同冻结核和同近似的 TMS：

```text
δ(¹H or ¹³C) = σ(TMS, same method) − σ(analyte, same method)
```

不允许：

- 用 PBE0 的 TMS 去定标 revDSD；
- 用溶液 TMS 定标一个未经说明的气相 CCSD(T) 值；
- 在看到每个实验峰后对每个原子单独加偏移；
- 隐藏参考物计算失败或方法不一致。

若以后建立经验线性标定，拟合参数只能来自明确的训练集，并必须在未参与拟合的化合物上验证。

## 8. J 耦合与多重性

J 主线候选为 `PBE/pcJ-3`，PBE0/pcJ-2 作为代表集基准。计算 J 后应使用完整或适当截断的自旋 Hamiltonian，而不只是给每个峰套一个经验 `s/d/t/m` 标签。

交换性 NH 的多重性容易被质子交换、氢键和弛豫抹平。原始计算 J 应保留供审计，但默认不强制把交换性 H 的 J 表现为清晰精细裂分。

## 9. 预期精度不是保证

完成所有门控、同方法参考和显式溶剂构象平均后，可将下列值作为努力目标：

- 芳香和其他非交换 `¹H`：峰数、顺序和多重峰形可望较接近；位移误差争取在零点几 ppm；
- `¹³C`：通常适合结构鉴定，但个别峰仍可相差数 ppm；
- NH/交换 H：可能相差数个 0.1 ppm，甚至超过 1 ppm；
- 峰宽、峰高、噪声、基线和溶剂峰：不属于静态 ORCA 的唯一第一性原理预测。

本包的晋级目标为非交换 `¹H` MAE ≤ 0.20 ppm、`¹³C` MAE ≤ 2.5 ppm。这是预先登记的验证阈值，不是当前结果，也不能保证对许多药物分子仍成立。药物常有更多构象、盐型、互变异构、聚集和慢/快交换，误差通常更难控制。

## 10. 仪器外观层

外观层可以基于实验元数据或预先登记的仪器模型加入：

- DMSO-d₆ 残余峰和水峰；
- Lorentzian/Gaussian 或 Voigt 线宽；
- 基线漂移、相位误差和数字滤波近似；
- 随机噪声；
- 积分归一化和动态范围。

这些参数必须独立保存、固定随机种子并标记为 `SYNTHETIC_ACQUISITION_APPEARANCE`。它们不能改变 clean calculated analyte peak positions，不能复制一条实验噪声轨迹，也不能用于制造“一模一样”的假实验谱。

## 11. 验证设计

- `F3_005` 苯甲酰胺：校准化合物，具有本地 Bruker 原始/处理数据 provenance；
- `F3_015` 中性 `(S)`-脯氨酰胺：留出效果检查，目前只有处理谱图目标和不完整采集元数据。

评价应同时报告：

- 非交换 `¹H` MAE、最大绝对误差、峰顺序；
- `¹³C` MAE、最大绝对误差和峰顺序；
- 交换 H 的单独误差；
- 峰数与多重性人工检查；
- revDSD 与 PBE0 方法差；
- 显式溶剂计量敏感性；
- 失败或未收敛任务比例。

只有完整计算、参考哈希、留出指标和具名复核人都存在时，才能将状态从 `DRAFT_NOT_VALIDATED` 晋级。当前尚未晋级。

## 12. 官方 ORCA 6.1 依据

- [NMR properties tutorial](https://www.faccts.de/docs/orca/6.1/tutorials/spec/NMR.html)
- [MP2 and DLPNO-MP2](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/MP2.html)
- [AutoCI coupled-cluster methods](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/autoci.html)
- [VPT2](https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/vpt2.html)
- [Thermochemistry and multi-temperature PrintThermoChem](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/thermochemistry.html)
- [Solvation models](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)
- [Parallel execution](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html)

---

## English rationale summary

There is no universally “highest” ORCA setting for solution NMR. Electronic-structure error is only one component; conformer populations, microsolvation, protonation, reference consistency, exchange, temperature, and acquisition processing can dominate. This package therefore uses revDSD-PBEP86-D4/2021/pcSseg-3/NoFrozenCore as a production candidate, PBE0/pcSseg-3 as a full-ensemble cross-check, and RI-MP2, TightPNO DLPNO-MP2, AUTOCI-CCSD(T), and VPT2 as restricted, manually gated diagnostics.

Same-method TMS referencing is mandatory. Boltzmann averaging is restricted to an identical explicit-solvent stoichiometry. Each nucleus uses a QRRHO free energy evaluated at its locked acquisition temperature from the same isotope-aware Hessian. VPT2 does not model exchange kinetics, and an analyte-only VPT2 correction is blocked from the main delta. The permitted composite is `delta_static_high + Delta_delta_vib_low`, where the low-level static/vibrational correction is computed for both analyte and TMS and all evidence is hash-bound. Gas-phase AutoCI is diagnostic only and requires its separate response outputs. Instrument-like noise, baseline, solvent peaks, intensity, and linewidth are generated only in a separate synthetic appearance layer. The pre-registered accuracy thresholds are validation targets, not achieved results or guarantees for drug-like molecules.

Solution-equivalence averaging is also pre-declared rather than inferred from whichever explicit-solvent cluster gives a convenient spectrum. For benzamide, the two ortho H/C sites and two meta H/C sites are averaged under the constitutional ring symmetry, with raw atom-resolved shifts and J values retained. The NH2 protons are not merged. Neutral `(S)`-prolinamide has an explicit negative declaration because its stereocentre can make CH2 protons diastereotopic.
