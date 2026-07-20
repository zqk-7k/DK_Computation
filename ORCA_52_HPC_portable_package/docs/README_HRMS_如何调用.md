# 52种化合物的理论HRMS数据：调用方法

## 结论与边界

这套流程为52种登记结构生成：

- 登记分子式的单同位素精确质量；
- 常见且保守的ESI完整前体/加合离子候选；
- 每个候选离子的精细同位素棒图数据；
- 以设定分辨率生成的高斯峰形CSV和PNG；
- 汇总CSV、制表符TXT、逐化合物TXT、验证JSON和SHA-256清单。

所有输出均标记为 **CALCULATED PREDICTION — NOT EXPERIMENTAL**。它们不是实测HRMS，不预测哪个离子一定出现，也不生成碎片离子、污染峰、背景噪声或仪器响应。论文中可把它们用于“calcd for …”的理论值；“found”必须来自真实质谱仪。

ORCA不参与这一部分。ORCA用于量子化学结构、NMR屏蔽和耦合计算；由分子式得到HRMS精确质量和天然同位素包络，应使用独立的质量/同位素卷积程序。

## 输入与结构核验

程序读取并交叉核验：

- `..\Figure3_52_ORCA_ordered_reference_states.sdf`
- `..\Figure3_52_ORCA_reference_state_registry.csv`

每个SDF结构必须与登记表中的Figure ID、分子式、净电荷和InChIKey完全一致，否则程序停止。

盐的处理规则：

- F3_016、F3_017：显式Cl⁻作为反离子排除，主候选为 `[M−Cl]+`。
- F3_021：两个显式Cl⁻作为反离子排除，主候选为 `[M−2Cl]2+`。
- 这是ESI去溶剂化/反离子处理模型，不是碰撞碎裂预测。
- 两性离子生成 `[M+H]+` 与 `[M−H]−`候选。
- 中性含N/O化合物通常生成 `[M+H]+`、`[M+Na]+`；具有明确可去质子化位点时另生成 `[M−H]−`。

当前52种结构包含H、C、N、O、F、P、S、Cl和I。程序还包含Br、Na和K的同位素表，以避免后续含溴结构无提示失败。Cl、Br和S均执行M/M+2比例自检；P和I按单一稳定天然同位素处理。

精确同位素质量和代表性天然丰度以NIST数据库为依据：

https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?all=all&ascii=ascii2&isotype=all

需要注意：某些元素的天然同位素组成随材料来源而变化，因此理论丰度不是任何具体样品的绝对实测丰度。

## 直接运行全部52种

在PowerShell中执行：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_hrms_formula.ps1"
```

若结果目录已经存在并确认要重新生成：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_hrms_formula.ps1" -Overwrite
```

只测试一个化合物，例如F3_006：

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "${PACKAGE_ROOT}\run_hrms_formula.ps1" -Only F3_006 -Overwrite
```

不要在正式全量结果目录中保留单化合物测试结果；测试后应使用不带`-Only`的全量命令重建。

## 输出位置

`results\HRMS_formula_predictions_CALCULATED_NOT_EXPERIMENTAL\`

主要文件：

- `Figure3_formula_HRMS_prediction_summary.csv`：全部候选离子和calculated m/z。
- `ALL_candidate_ions.txt`：同一汇总的TXT版本。
- `ALL_isotope_sticks.csv`及`.txt`：所有理论同位素棒图点。
- `isotope_mass_and_abundance_table.csv`：计算所用同位素质量和丰度。
- `VALIDATION_AND_PROVENANCE.json`：结构数量、输入哈希、软件版本和元素自检。
- `FILE_MANIFEST_SHA256.csv`：输入、程序和生成文件的SHA-256。
- `F3_xxx\F3_xxx_HRMS_summary.txt`：每个化合物的人类可读摘要。
- `F3_xxx\*_isotope_sticks.csv`：精细同位素棒图。
- `F3_xxx\*_gaussian_profile_R60000.csv`：理论连续峰形数据。
- `F3_xxx\*_gaussian_profile_R60000.png`：理论峰形图。

## 资源占用

这不是ORCA计算。运行器固定`OMP_NUM_THREADS`、`MKL_NUM_THREADS`、`OPENBLAS_NUM_THREADS`和`NUMEXPR_NUM_THREADS`为1，不会开启并行ORCA/MPI任务。运行时间通常为分钟级，内存和CPU负荷远低于NMR量子化学计算。

## 论文中的正确写法

允许：

> HRMS (ESI) m/z: calcd for C7H8NO [M+H]+ 122.0600; found [真实仪器结果].

不允许：

> HRMS found 122.0600.

如果没有真实仪器结果，不能把理论值写成“found”，也不能把生成的PNG作为实测原始谱图提交。
