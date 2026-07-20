# NMR 两层后处理与绘图模块

该模块只读取已经计算完成的标准 stick/transition CSV。它**不会启动 ORCA**，也不会修改 GOAT、DFT、shielding 或 coupling 作业。

## 科学边界

- Layer A：`CALCULATED MODEL - NOT EXPERIMENTAL`。只含输入的计算分析物跃迁，不加入溶剂、水、噪声、基线或杂质峰。
- Layer B：`SYNTHETIC APPEARANCE - NOT ACQUIRED`。从复数 FID 开始，加入有名称且有账本的残余溶剂/水组件、复高斯噪声、指数窗、零填充、FFT、残余相位和慢基线。
- 脚本永远不生成匿名“杂质峰”。所有数字表和图内均保留计算/合成标识，不能用作实测 NMR、纯度或身份鉴定证据。
- DMSO-d6 的 1H/13C 采集、噪声和基线参数从本地真实谱图审计中的经验分布确定性抽样。DMSO-d6 的 1H 残余溶剂峰和水峰也从本地高置信度记录抽样。
- D2O 目前没有本地校准，所有相关输出明确标为 `ASSUMED_NOT_LOCALLY_CALIBRATED`。
- 19F/31P 已保留直接核种接口，但没有本地对应核种采集校准时会使用明确标注的预设。

## 输入 CSV

必须包含以下列：

```text
transition_id,component_id,nucleus,shift_ppm,intensity,linewidth_hz,assignment,data_status
```

每条裂分线应是一行。`data_status` 只允许 `CALCULATED` 或 `CALCULATED_PREDICTION`。`transition_id`、`component_id` 和 `assignment` 不得为空。

测试夹具位于：

```text
tests/fixtures/nmr_standard_transitions_fixture.csv
```

## 单个谱图调用示例

```powershell
py -3 scripts\nmr_acquisition_appearance.py `
  --input results\NMR_predictions_CALCULATED_NOT_EXPERIMENTAL\F3_019\F3_019_1H_transitions.csv `
  --spectrum-id F3_019 `
  --display-name "compound F3_019" `
  --nucleus 1H `
  --solvent D2O
```

13C、19F 或 31P 只需改变 `--nucleus`。生产运行不要设置 `--acquired-points` 或 `--processed-points`；这两个参数只用于快速兼容性测试或有记录的人工覆盖。

默认输出到：

```text
results/NMR_acquisition_appearance_SYNTHETIC_NOT_ACQUIRED/<spectrum-id>/<nucleus>/
```

每个谱图目录包含：

- `layer_A_calculated_model.csv/.txt`
- `layer_B_synthetic_acquisition.csv/.txt`
- `layer_B_complex_FID.csv`
- `component_ledger.csv`
- `input_transitions_normalized.csv`
- `spectrum_metadata.json`
- `output_file_hashes.csv`
- `run_complete.json`
- `*_NMR_layers.svg`（可编辑主图）
- `*_NMR_layers.png`（600 dpi 预览/交付图）

## 测试

```powershell
py -3 -m unittest discover -s tests -p "test_nmr_acquisition_appearance.py" -v
```

测试覆盖确定性独立 seed、DMSO 命名组件、D2O 假设标记、13C/19F/31P 接口、PNG/SVG 生成和 SHA-256 清单。
