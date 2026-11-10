# W L3 补算能区与 FT-EXAFS 距离

目标是让 FEFF 原生输出覆盖 **10100–10300 eV**，然后导出 TXT 和图片。此目标由原生 `xmu.dat` 的第一列检验；扩大图片横坐标、外推或补零均不能替代补算。

## 输入设置

FEFF 的 `EGRID` 使用相对于其计算费米能级的能量，不接受将 10100、10300 直接当作该相对网格的端点。建议初始输入为：

```text
XANES 4.0 0.05 0.2
FMS 8.0 0
RPATH -1
EGRID
e_grid -160 130 0.5
POTENTIALS
...
```

`EGRID` 后面的网格行由下一个合法控制卡（这里为 `POTENTIALS`）结束。不要使用旧的 `EGRID 2 ...` 形式：当前固定版本的源代码指出该旧形式可能不工作。该新网格含 581 个水平能量点，并给复数能量积分路径留下余量；固定源码的 `nex=2000`，本包将用户网格上限设为 1800 点。

配置形式为：

```json
{
  "feff": {
    "egrid_by_edge": {
      "W_L3": [-160, 130, 0.5],
      "Co_K": [-160, 130, 0.5],
      "Fe_K": [-160, 130, 0.5]
    },
    "absolute_windows_eV": {"W_L3": [10100, 10300]}
  }
}
```

Co、Fe 也使用相同的相对能区，但按各自原生光子能量绘图。EXAFS 仍独立计算至 k=14 Å⁻¹，不把整个 EXAFS 能区强行交给 XANES 的 FMS 计算。

根据本地已有的上一轮原生 W L3 XANES，7 个含 W 结构的原始能区只有约 10201.56–10267.81 eV。这说明必须重新计算能量网格。新的几何结构可能改变阈值，所以上面的相对网格只是初始设置，最终以新输出为准。

## 自动验收及失败处理

`prepare_inputs.py` 在每个 XANES `job.json` 记录相对网格和绝对目标能区。`run_feff_job.py` 只在以下检查通过后标记任务完成：

1. 原生表格非空、数值有限、光子能量严格递增。
2. W L3 第一列最小值不高于 10100 eV，最大值不低于 10300 eV。
3. 原生间隔没有超过指定步长的异常缺口。

独立检查命令：

```bash
.venv/bin/python scripts/extend_feff_inputs.py --job jobs/PROFILE/MODEL/W_L3/site_NNNN/xanes
```

若不通过，保留原始尝试，在新 profile 中调整相对网格并重新运行。可由实际输出首点减去输入网格首点估计该次计算的相对/绝对能量偏移，再据此设置留有余量的端点。此估计只用于选择下一次求解网格，**不能用于把每条结果独立移到人为指定的边位**。同一谱线内的共同能量差也可用于交叉核对。若异常位移很大，应先检查势能计算、芯空穴、结构与软件日志。

FEFF 手册提示 FMS 在约 k>6 Å⁻¹ 时不再适合常规近边计算。因此默认正向端点为 +130 eV（约 k=5.84 Å⁻¹）；如调整需要超过约 +137 eV，本包要求明确审阅并启用 `allow_high_energy_fms`，不能仅为覆盖坐标而盲目扩大。远负能量端的输出由 FEFF 自身的复数能量与寿命展宽处理产生，不构造人工基线；它不能被解释为已经计算出全部束缚态激发。

与实验比较时若有参考样品，可采用有记录的统一能量校准。没有参考时保留 FEFF 原生能量，明确其绝对边位精度有限；不对齐每条白线峰来消除或制造状态差异。

## FT 图约 1.5 Å 的含义

未进行散射相移修正的 FT-EXAFS 幅值峰，不直接等于真实原子间距。实际 W–O 距离由结构坐标或有依据的散射路径拟合得到。W–O 的峰在 1.5 Å 左右可以作为待比较的谱学特征，但不能通过把 W–O 键设成 1.5 Å、平移 R 轴、删掉金属散射路径或调整参数来强制实现。

同一组比较采用相同的 k 区间、加权和窗口，先平均带符号的 χ(k)，再计算复数 FT。保存 R、幅值、实部、虚部以及所用窗口。W 第一壳全为 O 不代表完整 EXAFS 不应包含更远 Co/Fe 壳层或多重散射。

## 官方依据

- FEFF10 用户手册：XANES（第 70–71 页）、EGRID（第 76 页）与输出定义。https://times-software.github.io/feff10/assets/doc/feff10_users_guide.pdf
- 固定版本网格实现，允许用户定义负相对能量并处理复数能量路径。https://github.com/times-software/feff10/blob/d861cf8b64342e2d63a0aeb62c8b78faca7942c3/src/XSPH/phmesh2.f90
- 原生谱输出遍历完整水平能量网格。https://github.com/times-software/feff10/blob/d861cf8b64342e2d63a0aeb62c8b78faca7942c3/src/FF2X/ff2xmu.f90
- 网格点数上限。https://github.com/times-software/feff10/blob/d861cf8b64342e2d63a0aeb62c8b78faca7942c3/src/COMMON/m_dimsmod.f90
- Larch 的傅里叶变换定义及带散射相位的 χ(k) 表达式。https://xraypy.github.io/xraylarch/xafs_fourier.html
