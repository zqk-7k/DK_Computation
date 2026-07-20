# 可复制执行指令：Figure 3 的 52 个化合物 ORCA NMR 模拟

## 目标和科学边界

以 `Figure3_52_ORCA_reference_state_registry.csv` 和 `Figure3_52_ORCA_ordered_reference_states.sdf` 为唯一结构注册源，对 F3_001–F3_052 分别生成尽可能接近指定溶液 NMR 条件的计算谱。所有输出必须标为计算/合成结果，不得写成实测数据，也不得用于替代原始 FID、真实谱图或产物鉴定。

没有一条 ORCA 指令能直接生成与实验完全一致的整张谱。完整结果必须分成三层：

1. 分子层：构象、物种、屏蔽常数和 J 耦合；
2. 溶液层：隐式溶剂、显式氢键团簇、构象/物种加权；
3. 采集层：场强、线宽、峰形、溶剂/水背景、噪声、基线、相位和数字处理。

## 本地调用入口与强制串行策略

可执行入口不是本 Markdown，而是：

`${PACKAGE_ROOT}\run_serial.ps1`

完整调用命令见相邻项目中的 `README_如何调用.md`。默认是 dry-run，必须显式加入 `-Execute` 才允许进入计算，并且只有审核通过的 `condition_lock.json`、`stage_manifest.csv` 和 ORCA 输入齐全时才启动。

本机固定为一次一个化合物、一次一个阶段、一次一个构象。GOAT 使用 4 个 MPI 进程和 `NWorkers 4`，满足 ORCA 对 worker 数为 4 的倍数的要求；其他阶段默认 6 个 MPI 进程，OMP/MKL/OpenBLAS 各 1 线程。一般优化、频率和屏蔽使用 6×4,000 MB，J 耦合使用 6×4,500 MB，F3_006 的 X2C 或大型显式团簇使用 4×5,800 MB。并行资源数只影响耗时，不改变电子结构方法、构象覆盖或收敛阈值。

F3_019、F3_020、F3_021 使用 D2O 主条件，DMSO-d6 只作为可选第二条件；F3_035、F3_043、F3_044 使用 D2O。两种溶剂中的物种和构象总体不得混合平均。

## 0. 运行前硬性检查

对每个 F3 编号建立不可变的 `condition_lock.json`，至少记录：

- figure_id、名称、Isomeric SMILES、InChIKey、结构哈希；
- 电荷、多重度、手性、质子化/互变异构状态、反离子；
- 实际氘代溶剂、样品浓度、温度、pH/pD；
- 仪器 ^1H 标称频率、核种、脉冲程序、扫描数、谱宽、采样点、采集时间、d1；
- 参考物和参考方式；
- 原始实验记录的来源。

若这些信息来自本方案的建议表而不是实验记录，写 `condition_status = PROPOSED_NOT_VERIFIED`。不得把建议条件改写成“实际测试条件”。

## 1. 物种和构象

### 1.1 普通中性分子

1. 用多个独立准备并哈希锁定的起始三维构象/显式溶剂排布；ORCA 6.1 的 `RandomSeed false` 用于确定性运行，但它不是一个可编号的随机种子，因此不把三次搜索写成“seed 1/2/3”。
2. 使用 ORCA `GOAT-ENTROPY XTB ALPB(<solvent>)` 搜索构象；柔性分子至少从 3 个不同起始几何逐个进行 GOAT 搜索。
3. 合并、按重原子 RMSD 和能量去重；GOAT 搜索窗建议 12 kcal mol^-1。
4. 对进入 DFT 的构象使用 `r2SCAN-3c/CPCM(<solvent>) Opt TightOpt Freq`。
5. 仅保留无虚频的极小值。GOAT 使用锁定的 `ensemble_search_temperature_K`；最终 QRRHO/Boltzmann 权重则分别按 ^1H、^13C、^19F 或 ^31P 的实际/明确假设采集温度重算。至少保留所有 DFT 相对自由能不高于 3.0 kcal mol^-1 的构象，并用 5.0 kcal mol^-1 窗口做敏感性分析。

### 1.2 盐、两性离子和 pH 依赖物种

- F3_016、F3_017：分别比较阳离子连续介质模型、氯离子接触离子对、DMSO 分隔离子对。
- F3_021：比较双阳离子连续介质模型、两个氯离子的接触/溶剂共享模型；胍基互变异构和质子化必须锁定。
- F3_019、F3_020、F3_023：由实际溶剂和 pH/pD 决定中性或两性离子模型。尤其 F3_020 不能因为图中画了中性形式就默认其在 D2O 中仍为中性。
- F3_035、F3_043、F3_044：酸碱状态对水溶液谱的影响可能大于泛函误差，必须建立 pD 对应的物种方案。

不同质子化状态、不同反离子状态以及含不同数目显式溶剂分子的团簇，不能仅凭彼此的 ORCA 自由能直接做 Boltzmann 混合；它们具有不同化学组成和溶剂化学势。没有实验物种比例时，应分别输出为情景谱和不确定性包络。

### 1.3 与主溶剂一致的显式溶剂模型

外层连续介质和显式第一溶剂化层都必须与锁定的主溶剂一致：

- **DMSO-d6 主条件**：外层用 DMSO；对 NH/OH 位点建立 1–2 个显式 DMSO 的代表性氢键模型。痕量水只作为另一个固定水活度/化学计量的敏感性情景，不可无权重地强制并入主模型。
- **D2O 主条件**：外层电子连续介质用 `WATER`，第一层使用显式水合团簇，不加入 DMSO。F3_019/020 依次测试 6、8、10 个水，F3_021 测试 8、10、12 个水，F3_035/043/044 测试 4、6、8 个水；以加权位移变化小于 ^1H 0.05 ppm、^13C 0.5 ppm、^31P 1.0 ppm 作为收敛目标。
- 每一种固定团簇化学计量至少准备 10 个起始氢键构型，并逐一运行；不同显式溶剂计量、不同质子化状态和不同反离子组成绝不直接做 Boltzmann 混合。

显式溶剂原子必须在 atom map 中标为 `environment_atom = true`，并从最终溶质 NMR 信号中排除。普通 H2O 或普通 DMSO 的计算氢绝不能直接充当实验中的 HOD/水峰或 2.50 ppm 背景峰；这些实验背景只能在独立、标记清楚的外观后处理层中加入。

在 Born–Oppenheimer 电子结构层面，H 和 D 的核电荷相同；但 D2O 主条件的频率/QRRHO 需要按 ORCA 6.1 坐标语法为环境水及发生交换的位点指定 D 质量（约 2.014 u），并保存 isotope map 与 SHA-256。若兼容性测试不能确认该表示，只能另标为 `H2O-mass approximation` 并做敏感性分析，不能冒充 D2O 同位素自由能。无论采用哪种表示，都不得把显式溶剂或快速交换位点的氢积分计入产物积分。

## 2. 几何优化和频率

每个构象按主条件选模板：DMSO-d6 使用 `ORCA_input_templates/01_opt_freq_DMSO.inp`，D2O 使用 `ORCA_input_templates/01_opt_freq_D2O_primary.inp`。把电荷、多重度、并行数、内存、温度和坐标文件替换为该任务的锁定值。

验收条件：

- `ORCA TERMINATED NORMALLY`；
- 最终结构无虚频；极软扭转模需要人工审核而不是静默删除；
- 没有断键、质子转移到错误位点、盐的非物理解离或显式溶剂逃逸；
- SCF 和几何收敛阈值满足输入设置。

## 3. 屏蔽常数

对所有保留构象使用 `PBE0/pcSseg-2/GIAO/CPCM(<solvent>)`。DMSO-d6 使用 `ORCA_input_templates/02_shielding_H_C_F_P_explicit_indices.inp`，D2O 使用 `ORCA_input_templates/02_shielding_D2O_primary_explicit_indices.inp`。至少计算 ^1H 和所有存在的 ^13C；F3_003、F3_037、F3_039 另算 ^19F；F3_043–F3_052 另算 ^31P。

计算的原始屏蔽标尺与实验定标必须分开记录。^1H/^13C 可以保留同层级 TMS 计算得到的 raw scale，但最终实验尺度必须匹配实际测试的 reference method：DMSO-d6 实谱可按残余 DMSO（约 2.50/39.52 ppm）校正；D2O 必须按实际 DSS/TSP、间接或外标方法校正，不能默认使用不溶于水的 TMS。^19F 和 ^31P 各用独立、经审核的参考体系，绝不套用 TMS。

每一种 `核种 + 方法 + 溶剂域 + 参考方式` 只建立一次 reference job，并在 `reference_job_registry.json` 中保存输入、输出和方法 SHA-256。目标分子的 condition lock 只能引用已审核的 reference ID；不得逐分子改参考值以追近实验峰位。

同时保留两套位移：

1. `raw_reference_shift`：直接由同层级参考屏蔽常数换算；
2. `calibrated_shift`：只允许使用预先锁定、与当前计算层级和溶剂域匹配的外部标定方程。

禁止逐个分子或逐个峰拟合到文献/实测值。任何标定方程必须先于目标集锁定，并报告训练集、适用域、RMSE 和残差。

直接参考换算同时保存严格式和常用近似式：

\[
\delta_i=\frac{\sigma_{ref}-\sigma_i}{1-\sigma_{ref}\times10^{-6}}
\quad\text{和}\quad
\delta_i\approx\sigma_{ref}-\sigma_i
\]

在全批计算前设置方法精度门。代表集至少包括 F3_003（^19F）、F3_006（I）、F3_016/F3_021（盐）、F3_022 或 F3_023（柔性/两性离子）、F3_030 或 F3_033（慢酰胺旋转）及 F3_046 或 F3_050（^31P）。比较 PBE0/pcSseg-2 与至少一个更大基组方案（如 PBE0 或 TPSSh/pcSseg-3）；小型代表体系可再做更高层级敏感性计算。方法必须依据预先定义的验证指标统一锁定，不能为每个化合物选择“看起来最像”的方法。

### F3_006 的碘

`pcSseg-2` 只覆盖到 Kr。F3_006 不得给 I 直接套用不存在的 pcSseg-2 参数。非相对论主分支使用 `02a_F3_006_shielding_mixed_basis_main.inp`：H/C/N/O 用 pcSseg-2，I 显式使用 def2-TZVP/def2-ECP；相对论敏感性分支使用 `02b_F3_006_X2C_sensitivity.inp`。J 耦合相应使用 `03b_F3_006_J_iodine_special.inp` 和 `03c_F3_006_J_X2C_sensitivity.inp`。四个模板都必须先通过 ORCA 6.1 小型语法/基组/性质兼容门；不同方法分支各用自己的同层级 reference ID，不能共用一个参考计算。两套差异进入不确定性报告。

## 4. J 耦合

使用 `PBE/pcJ-3/CPCM(<solvent>)` 作为建议的高质量方案；DMSO-d6 使用 `ORCA_input_templates/03_J_coupling_explicit_indices.inp`，D2O 使用 `ORCA_input_templates/03_J_coupling_D2O_primary_explicit_indices.inp`。现有本地流程已经验证 PBE0/pcJ-2 可运行，因此先在代表集上盲法比较 PBE/pcJ-3 与 PBE0/pcJ-2，再为全批统一锁定一种；不得按单个谱的外观临时切换。至少计算所有 H–H 耦合；含 F 的分子增加 H–F、C–F、F–F；含 P 的分子增加 H–P、C–P、P–P。F3_006 的 I 不在 pcJ-3 覆盖范围内，必须使用 `03b_F3_006_J_iodine_special.inp` 的显式 I 基组/ECP，并先通过小型语法/基组兼容门，再做独立 X2C 敏感性计算。

对超大自旋系统，不得试图一次对全部质子做不可行的 2^N 精确对角化。按耦合图和强耦合判据切分局部自旋系统：

- 小型、强耦合子系统：精确自旋 Hamiltonian；
- 弱耦合远程网络：一阶分裂；
- 对两种处理边界做敏感性比较。

^13C{^1H} 常规去偶谱不显示 C–H 分裂；^19F 和 ^31P 应分别输出去偶与未去偶版本，以匹配真实采集程序。

## 5. 构象平均和动态等价

对相同化学组成的构象，用锁定的 QRRHO 权重分别平均屏蔽常数和 J。保留每个单构象值、相对自由能和权重。

不要仅根据二维结构自动写 `NMREquiv`。只有在原子映射、分子对称性和溶液中的快速动力学共同支持时才可设为等价；慢酰胺旋转、环翻转、苄位非对映异位和盐的离子配对应保留可能的非等价或次要转子信号。

## 6. 从 ORCA 参数生成谱

^1H 使用 400.13 MHz；对应的建议显示频率为 ^13C 约 100.61 MHz、^19F 约 376.5 MHz、^31P 约 162.0 MHz。ORCA 的 `orca_nmrspectrum` 可用于小自旋系统检查；主生产谱建议由可审计的 Python 后处理程序从逐核位移和 J 矩阵生成。

保存：

- 逐原子屏蔽常数、原始/标定位移；
- 完整 J 矩阵；
- 构象自由能和权重；
- 理想 stick spectrum；
- 无噪声展宽谱；
- 实验外观模拟谱；
- 每条连续曲线的 TXT/CSV；
- 完整参数 JSON、后处理噪声随机种子和 SHA-256（此处的随机种子只指合成外观层，不宣称为 GOAT 的随机种子）。

`orca_nmrspectrum` 建议从 `NMRCoal = 0.10 Hz` 开始，避免人为合并精细线；最终视觉合并由真实线宽和数字分辨率完成。`NMREquiv` 的原子编号必须来自当前 XYZ/atom map，不能手填猜测。

## 7. 实验外观后处理

使用 `experimental_appearance_postprocessing.json`。优先在时域生成复数 FID：

1. 由化学位移/J 转换成各跃迁频率和强度；
2. 用 T2* 对应的 Lorentzian 衰减和小量 Gaussian 场不均匀性生成 Voigt 线形；
3. 加入有固定随机种子的复高斯热噪声及较弱相关噪声；
4. 加入锁定的全局频率漂移、极小相位误差和低阶基线；
5. 按规定窗函数、零填充、Fourier 变换和自动相位/基线流程处理。

必须同时输出：

- `CALC_ONLY`：不含噪声、基线、溶剂、水或杂质；
- `SIMULATED_ACQUISITION`：含全部已声明后处理成分；
- `COMPONENTS`：产物、溶剂、水、允许杂质、白噪声、相关噪声、基线各自单独的曲线。

### DMSO-d6 背景

- ^1H 残余 DMSO：中心 2.500 ppm；高分辨条件下可建模为 DMSO-d5 的细小 quintet（约 1.9 Hz 的 H–D 耦合），若实际分辨率不足则呈近似 singlet。
- DMSO 中水：3.33 ppm 只是 298 K 附近的名义值。峰位、宽度和强度随水量、温度、浓度及交换而变；没有实测空白时用不确定性分布，不能固定成每个样品完全一样。
- ^13C DMSO-d6：约 39.52 ppm 的 C–D septet，约 21 Hz；相对线强采用 1:3:6:7:6:3:1，并根据数字分辨率和线宽决定是否可见精细分裂。

### 噪声、基线和额外小峰

- 噪声必须在 FID 层加入，并以“最大非溶剂产物峰”为基准定义目标 SNR，不能对每个 ppm 点直接添加肉眼好看的独立抖动。
- 允许白噪声加弱相关噪声；每个谱使用由 figure_id 派生、可重复的随机种子。
- 峰位波动分为全谱共同漂移和化学环境波动。仪器漂移应是全局量；水和交换性质子可有更大的独立波动；不得任意移动稳定的非交换产物峰。
- 额外杂质峰只能来自真实实验记录、溶剂空白或预先声明的候选库。没有证据时，不生成“未知杂质”。若仅为算法压力测试，可另生成 `UNCERTAINTY_DEMO_ONLY`，且不得进入表征文件。
- ^13C{^1H} 峰高不能设为所有碳等高，也不能解释成碳数。应根据脉冲角、d1、T1、NOE 和碳类型建立经验响应范围；没有 T1 数据时只作为显示模型。

## 8. 条件特例

- D2O：不加入 2.50/3.33 ppm；加入约 4.79 ppm 的 HOD 背景，并从最终 ^1H 谱中删除快速交换的 NH/OH。使用 CPCM(water) 和与 pD 对应的物种。
- CD3OD：不加入 DMSO 背景；使用约 3.31 ppm 的残余 CD3OD 和约 4.87 ppm 的水/HOD背景，并考虑交换性质子消失或显著变宽。
- CDCl3：不加入 DMSO 背景；使用 7.26 ppm 的 CHCl3 残余峰，水峰位置取决于样品条件。
- F3_035、F3_043、F3_044 不含碳，不生成虚构的 ^13C 谱。

## 9. 质量控制与停止条件

每个化合物必须通过：

- 结构、手性、电荷、反离子和原子映射一致；
- 所有采用的 DFT 构象正常结束且为极小值；
- 构象权重和信号积分守恒；
- 对称/动态等价有依据；
- 溶剂/水环境原子未混入产物积分；
- DMSO/D2O/CD3OD/CDCl3 背景与条件一致；
- 随机种子、后处理参数和组件曲线齐全；
- 图和数据永久标为 `SIMULATED / NOT ACQUIRED`。

若实际溶剂、物种、pH/pD 或参考方式不明确，应输出情景比较并标记 `BLOCKED_FOR_EXPERIMENTAL_METADATA`，不得悄悄选择一个最像文献的结果。

## 10. 最佳校准路径

已只读盘点本地 `处理过的谱图`。详细结果在 `real_spectrum_calibration/processed_spectra_inventory.csv` 和 `benzamide_calibration_notes.md`。当前没有可独立验证的原始 Bruker vendor FID；`.mnova`、PDF、TIFF 和处理后 XY 不得直接改称 FID。

29 条嵌入记录给出了一套稳定的实验室采集指纹：16 条 ^1H 均为 NS=8、D1=1 s、LB=0.3 Hz、SI=65536；13 条 ^13C 均为 NS=512、D1=2 s、LB=1.0 Hz、SI=32768；场强约 400.13/100.61 MHz。苯甲酰胺真实处理谱进一步给出 ^1H 非交换线宽约 3.3 Hz、交换峰约 18 Hz，以及 ^13C 产品峰 1.77–4.45 Hz。`experimental_appearance_postprocessing.json` 已据此更新。

如目标是“最接近本实验室的真实测试”，最有效的信息不是再增加随机噪声，而是继续提供同一台仪器、同一探头和同一采集模板下的：

1. 一个 DMSO-d6 空白 FID；
2. 苯甲酰胺原始 Bruker 实验目录，或由 Mnova API 成功导出且可校验的复数 FID；
3. 一个离子型样品、一个含 P 样品的真实 FID；
4. 对应的采集参数和处理历史。

利用这些数据只校准仪器外观层（噪声功率谱、线形、相位、基线、DMSO/水强度分布），不要逐峰校准 52 个目标分子的化学位移。锁定参数后再盲法生成全套结果。

## 11. 主要方法依据

- ORCA 6.1 NMR 教程：https://www.faccts.de/docs/orca/6.1/tutorials/spec/NMR.html
- ORCA 6.1 NMR 手册：https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/nmr.html
- ORCA 6.1 隐式溶剂说明：https://www.faccts.de/docs/orca/6.1/tutorials/prop/cpcm.html
- ORCA 相对论方法说明：https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/relativistic.html
- 常见残余溶剂/杂质 NMR 基准：Fulmer et al., Organometallics 2010, DOI 10.1021/om100106e；Gottlieb et al., J. Org. Chem. 1997, DOI 10.1021/jo971176v。
