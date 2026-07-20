# ORCA NMR 参考标准工作流（当前仅 DRAFT）

## 当前状态

本层只准备“与目标 shielding 方法匹配的计算参考模型”。截至 2026-07-20：

- 6 个参考模型、输入、原子映射、声明、清单和门控已经生成；
- 默认命令只做 dry-run；
- `reference_model_protocol_DRAFT.json` 仍为 `DRAFT_NOT_APPROVED`；
- 没有运行任何参考标准 ORCA 作业；
- 计算 raw shielding、实验 reference convention 和谱图后处理仍是三个独立概念。

这些计算结果以后也只能标为 **CALCULATED MODEL – NOT EXPERIMENTAL DATA**。

## 已准备的模型

| 溶剂域 | 核种 | 计算参考模型 | 重要限制 |
|---|---|---|---|
| DMSO-d6 | ^1H、^13C | TMS / CPCM(DMSO) | 这是计算绝对 shielding 的零点模型，不表示实验样品加入了 TMS；实际谱图也可能以残余 DMSO 定标。 |
| DMSO-d6 | ^19F | CFCl3 / CPCM(DMSO) | 必须另行记录实验采用内标、外标还是间接 CFCl3 标度。 |
| DMSO-d6 | ^31P | H3PO4 + 2 DMSO 微溶剂化模型 | **不等同于 85 wt% H3PO4 外标**。 |
| D2O | ^31P | H3PO4 + 4 D2O 微溶剂化模型 | **不等同于 85 wt% H3PO4 外标**，也不能重现浓酸物种分布和体磁化率。 |
| D2O | ^1H、^13C | DSS^−/Na^+ + 6 D2O 候选模型 | 仅在真实测试确实采用匹配 DSS 条件时使用。 |
| D2O | ^1H、^13C | TSP^−/Na^+ + 6 D2O 候选模型 | TSP 与 TSP-d4 不可默认为同一结构。 |

DSS、TSP 和仪器/IUPAC 兼容的间接参考是三条互斥路线。批准门控时必须选定一条，并写明真实标准品、浓度、pD、温度、内/外标方式或间接参考参数。

## 已实施的安全校验

- TMS 的参考位点为 12 × ^1H、4 × ^13C；DSS/TSP 均严格为三甲基硅基的 9 × ^1H、3 × ^13C，排除了 Si–CH2 链段。
- `Nuclei` 使用 ORCA 的 1-based 索引；环境 DMSO/D2O 和 Na^+ 不会被请求为参考核。
- D2O 及 DMSO-d6 环境中的 D 质量以 `M=2.01410177812` 写入频率输入。
- 多组分起始簇经过跨组分碰撞检查；这只说明“没有明显原子重叠”，不说明已经找到正确的溶剂化构象。
- opt/freq 使用 `r2SCAN-3c/CPCM`；shielding 使用 `PBE0/pcSseg-2/AutoAux/RIJCOSX/GIAO/CPCM`，方法指纹和实际输入均锁定 SHA-256。
- shielding 的输入文本与最终 `optimized.xyz` 共同形成 bundle hash；输出、可执行文件及每个产物也会记录 SHA-256。
- opt/freq 必须正常结束、为零虚频、ORCA 报告 298.15 K 和 100 cm^-1 Grimme QRRHO 后，才允许进入 shielding。
- 一次只运行一个模型、一个阶段；MPI 为 6 ranks，超过持续 CPU 阈值时完成当前阶段后停止，不启动下一个阶段。

## 1. 重新生成 DRAFT（不会运行 ORCA）

仅当 `reference_jobs` 中没有任何 `.out`/`.gbw` 计算产物时，才可重建：

```powershell
cd ${PACKAGE_ROOT}
py -3.13 scripts\prepare_reference_jobs.py --project . --replace-draft
```

脚本会更新：

- `reference_jobs/<model_id>/start.xyz`
- `reference_jobs/<model_id>/reference_atom_map.csv`
- `reference_jobs/<model_id>/opt_freq.inp`
- `reference_jobs/<model_id>/shielding.inp.template`
- `reference_jobs/<model_id>/model_declaration.json`
- `data/reference_job_manifest_DRAFT.csv`
- `jobs/reference_model_protocol_DRAFT.json`
- `config/reference_job_registry.json`

## 2. 只读 dry-run

```powershell
py -3.13 scripts\run_reference_serial.py --project .
```

预期输出为：

```text
DRY RUN: validated 6 reference models; gate=DRAFT_NOT_APPROVED; ORCA was not started.
```

可只检查一个模型：

```powershell
py -3.13 scripts\run_reference_serial.py --project . --only REF_DMSO_TMS
```

## 3. 审核和批准要求

当前不要执行 `--execute`。总门控的一个“批准”不能替代以下逐项审核：

1. 锁定实验真实参考方式、标准品/同位素体、pD、温度、浓度和内标/外标/间接参考参数。
2. 检查每个模型的结构、质子化、原子索引和初始簇；显式簇最好再经过独立的 xTB/CREST/GOAT 构象检查。
3. D2O 的 ^1H/^13C 只能选 DSS、TSP 或间接参考中的一条。若选择间接参考，DSS/TSP 两个 ORCA 模型均应拒绝。
4. 每个分子模型必须标为 `APPROVED_FOR_REFERENCE_CALCULATION` 或 `REJECTED_REFERENCE_MODEL`；拒绝项必须写明原因。
5. 选中的模型声明、manifest 行和 gate decision 状态必须一致；更新声明后需重算其 SHA-256，再更新 manifest hash 和 gate hash。
6. 总门控的三个全局检查、每个获批模型的四个审核字段、路线决定、审核人和时间均需完整。

工作流有意不提供“自动批准”命令，以避免只改一个状态字段便绕过哈希和科学审核。

## 4. 最终获批后才可串行执行

```powershell
py -3.13 scripts\run_reference_serial.py `
  --project . `
  --gate jobs\reference_model_protocol_APPROVED.json `
  --execute
```

运行器会在 `reference_runs/<model_id>/<attempt_id>/` 保存独立阶段，在 `logs/reference_serial/jobs/` 保存状态。任何 DRAFT、状态不一致、路径越界、输入/模型/方法哈希变化、外部 ORCA/MPI/xTB 进程、资源门控失败或虚频，都会在启动下一阶段前中止。

## 5. 计算完成后的处理规则

- 保留每个获批参考原子的 raw isotropic shielding，不可只留下一个手填数字。
- 对声明为等价的参考位点，按核种报告 `n、算术平均值、总体标准差、最小值、最大值`。
- 不得跨模型、同位素体、溶剂域、构象或实验 processing convention 混合平均。
- 最终 `reference_job_registry.json` 必须绑定实际 shielding 输入、`optimized.xyz`、bundle、ORCA 输出及 ORCA 可执行文件的 SHA-256。
- 只有随后独立批准的 registry 才能被 52 个化合物的 property builder 使用。

## 尚未解决的科学事项

- D2O 中 DSS/TSP/间接参考的真实路线尚未从原始 acquisition/processing 记录确认。
- 两个 ^31P 微溶剂化 H3PO4 模型仅是声明清楚的计算近似，不能代替真实 85 wt% H3PO4 外标。
- 含碘/X2C 分支的方法与普通 PBE0/pcSseg-2 分支不同；在 F3_006 的 X2C gate 和匹配参考模型获批前，property builder 应继续 fail closed。
- 碰撞筛查不等于构象收敛。若这些显式参考簇要承担定量校准，应先比较多个物理合理构象，而不是依赖单一起始排布。
