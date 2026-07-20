# Figure 3 五个结构的最终参考态

这些结构用于消除二维结构歧义并准备 ORCA 输入。这里的“最常见/较稳定参考态”是有条件的：盐型、质子化状态和两性离子比例会随溶剂、pH、浓度及固态形式变化。因此，下列选择不能替代对实际分离样品的确认。

| 编号 | 采用的明确结构 | 关键说明 |
|---|---|---|
| F3_016 | L-serinamide hydrochloride，S | 单盐酸盐；显式 `NH3+` 和 `Cl−`；与作者 SI 的盐酸盐名称一致。 |
| F3_017 | L-leucinamide hydrochloride，S | 单盐酸盐；显式 `NH3+` 和 `Cl−`；与作者 SI 的盐酸盐名称一致。 |
| F3_019 | L-asparagine zwitterion，S | 采用近中性水相常见的 `NH3+`/`COO−` 两性离子；原图本身没有证明 L/D。 |
| F3_021 | L-argininamide dihydrochloride，S | α-铵盐和胍鎓形成二价阳离子，配两个 `Cl−`；是常见分离盐型，但源稿未证明实际样品盐型。 |
| F3_023 | 3-(2-hydroxybenzamido)-L-tyrosine zwitterion，S | 保留源图的水杨酰胺连接；结构中没有硝基，因此不再称为 nitrotyrosine。该精确身份和手性仍需真实样品证据确认。 |

## 文件用途

- `Figure3_full_reference_states.*`：完整 Figure 3 的 PNG、SVG 和 ChemDraw CDX。
- `individual_structures/`：五个结构各自的 PNG、SVG、CDX、CDXML 和 MOL。
- `Figure3_52_ORCA_ordered_reference_states.sdf`：按 F3_001–F3_052 排序、带元数据的 ORCA 建模结构源。
- `Figure3_52_ORCA_reference_state_registry.*`：分子式、SMILES、InChI、InChIKey、手性和建模备注。
- `F3_016_017_019_021_023_reference_state_manifest.*`：本次五个结构的决策清单。

ORCA 计算前仍需按实验条件选择溶剂模型、离子对初始构型、显式水、构象搜索策略及相对布居；不能把二维参考态直接等同于实验 NMR 中的唯一微观状态。
