# 补算结果的 TXT 和图片

`scripts/export_supplement.py` 读取实际 FEFF 结果经 Larch 处理得到的 `COMPLETE_*` 目录。只有全部位点计算完整、校验值一致时才导出正式图。它不会生成缺失的曲线，也不会使峰位、强度或不同元素的响应符合预设结论。

运行方式：

```bash
.venv/bin/python scripts/export_supplement.py --config config.json --result-dir /path/to/COMPLETE_result
```

输出位于该结果目录中的 `FINAL_EXPORT_时间戳/`。原结果保持不变。

- `txt/means/`：每个模型、每条吸收边的全位点平均数据。包括 XANES、带符号 χ(k)、k¹/k²/k³χ(k)、FT 的实部/虚部/幅值、WT 的实部/虚部/幅值。
- `txt/sites/`：每个独立吸收原子的原生/公共网格数据及变换数据，与平均数据分开。
- `txt/means/.../WT_Cauchy_long.txt`：方便 Origin 等软件导入的五列 WT 数据：k、R、实部、虚部、幅值；同时保留矩阵格式 TXT。
- `txt/comparisons/`：每组模型的三曲线合并表，可直接复绘。W XANES 必须由真实计算覆盖 10100–10300 eV；缺少范围时停止导出并要求重新计算，不能补零、拼接虚构背景或外推。
- `figures/`：每组比较均生成三行 W/Co/Fe 的 XANES、FT 和 k³χ(k)，以及三行三列 WT。全部提供 PNG、PDF、SVG 和 TIFF；默认图片 600 dpi，PDF/SVG 保留可编辑文字。Arial 可用时使用 Arial；否则使用 DejaVu Sans 并在报告中记录。
- `export_report.json`：来源、方法、字体、对比指标、FT 峰区网格最大值和解释限制。
- `export_inventory.json`：输出文件的大小和 SHA-256。

TXT 使用 UTF-8、制表符分列、`#` 开头的注释。最后一行注释为列名，包含单位。它保留数值，不添加噪声，不对每条曲线单独归一化，不把峰对齐。谱图的纵轴为 FEFF 的计算吸收强度；不同吸收边的绝对强度不用于直接比较活性。FT 使用 k² 权重的幅值作为主图，并保留 k¹/k²/k³ 的完整数字表。

## 三组工况与模型标签

当前构型是对实验条件的候选解释，默认图例为 `Hydroxylated`、`Deprotonated (no EG)` 和 `EG bound at W/Co`。静态几何结构本身没有建立 OCP 或 1.6 V 的对应关系，因此图例不会直接标注这些电位。实验条件可在说明文字中列为待解释的目标；确切电位对应需要额外的界面热力学/恒电位依据。

`config.json` 的配置形式为：

```json
{
  "export": {
    "w_energy_window_eV": [10100, 10300],
    "dpi": 600,
    "comparisons": [{
      "id": "comparison_identifier",
      "title": "Candidate structural-state comparison",
      "series": [
        {"model": "exact_model_id_1", "label": "Hydroxylated"},
        {"model": "exact_model_id_2", "label": "Deprotonated (no EG)"},
        {"model": "exact_model_id_3", "label": "EG bound at W"}
      ]
    }]
  }
}
```

三个模型的 W、Co 和 Fe 结果都必须存在；缺少某个状态时不会复制另一个模型代替。

## 峰位和差谱

FT 图横轴为 **未作相位校正的 R**。约 1.5 Å 的峰可以对应更长的实际 W–O 键，不能把结构中的键长强行设成 1.5 Å。导出器仅报告 1.0–2.5 Å 区间内幅值最大的原网格点；落在区间边界时会记录提醒。该数字不是键长拟合结果，也不包含配位数或置信区间。真实键长由松弛结构的配位检查报告给出。

XANES 差异指标来自相同能量窗口上的原始平均 μ，作为比较量，不换算成氧化态或活性位点占比。WT 由同一组 EXAFS 计算，不能视为独立实验验证。

默认不作能量平移。若有独立参考校准，可设置 `energy_reference_shift_eV`，但仅允许同一吸收边共用一个值，并必须同时提供 `energy_reference_evidence` 的 `description` 和可读取的 `file`。不允许对各状态分别对白线峰进行平移。

## 程序验证边界

`tests/test_export_supplement.py` 在系统临时目录生成明确标注 TEST ONLY 的解析测试数据，检查 TXT 数值保真、四种图片格式、缺失范围拒绝、源数据篡改拒绝、缺失状态拒绝及无依据电位标签拒绝。这些程序测试不是 FEFF 计算，不进入正式结果包。
