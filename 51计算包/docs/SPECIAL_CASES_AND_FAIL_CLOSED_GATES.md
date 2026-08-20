# 特殊化合物与 fail-closed 路由

“Fail-closed”表示：在关键科学设置未确认前，流程应停止并生成明确的阻断记录，而不是用默认模板悄悄继续。

## F3_006：2-碘苯甲酰胺

碘的标量相对论、基组和 NMR picture-change 处理会影响屏蔽常数。普通轻元素 `pcSseg-3/pcJ` 模板不能直接视为可用。生产 property 作业前必须：

1. 明确选择并固定 X2C、ZORA 或 ECP 路线；
2. 检查 I 及所有元素的 orbital/NMR/J/辅助基组覆盖；
3. 检查 GIAO 和 picture-change 设置是否与所选相对论 Hamiltonian 兼容；
4. 使用 ORCA 6.1.1 运行最小输入和参考物 pilot；
5. 保存输入、基组展开、版本和输出哈希；
6. 由方法负责人批准后才能生成生产 shielding/J。

在门控通过前，可以进行明确标记的几何/采样试运行，但不得输出与普通轻元素分支同等级的最终移位。

## F、Cl 和 P：基组预检

- F：`F3_003`, `F3_037`, `F3_039`；
- Cl：`F3_016`, `F3_017`, `F3_021`；
- P：`F3_043`–`F3_052`。

对每个 stage 检查 orbital basis、RI/JK auxiliary basis、NMR basis 和 J basis。ORCA 能读入关键词并不等于所有元素实际使用了预期基组；日志中的 basis summary 必须解析并入 gate receipt。含 F/P 的体系可选输出 `19F/31P`，但只有在核种、参考和实验需求明确时才启用。

## 无碳化合物

`F3_035`（sulfamic acid）、`F3_043`（phosphoramide）和 `F3_044`（phosphorodiamidic acid）不含碳。其 `13C` 状态必须为 `NOT_APPLICABLE`，不能生成空图后写成“无峰计算结果”，更不能伪造 13C 谱。可考虑 `15N`、`31P` 或正交表征，但应成为单独核种任务。

## D2O 与交换性氢

D2O 路由涉及 `F3_019`, `F3_020`, `F3_021`, `F3_035`, `F3_043`, `F3_044`。静态计算可给出某一微观结构中 NH/OH 的 shielding，但实验中这些氢可能发生 H/D 交换、严重展宽或完全消失。因此：

- 把交换性 H 单独标注；
- 不将静态 NH/OH 峰存在性作为普通 `1H` 谱的强制预测；
- J 表可保留相关值，spin simulation 默认不加入其相干裂分；
- 显式水数目不能替代 pD、交换速率或水活度信息。

## F3_016/F3_017：盐与 cation-only

L-serinamide hydrochloride 和 L-leucinamide hydrochloride 的 full-salt（含 Cl⁻）与 cation-only 是不同组成。二者必须具有不同 species hash、场景、参考和标签，不能按普通构象混合。默认生产物种由样品实际盐型决定；未确认时结果标为 provisional。

## F3_019/F3_020/F3_021：物种确认

- `F3_019`：确认 L 构型以及实验条件下的两性离子状态；
- `F3_020`：registry 的中性 glutamine 图与 D2O 工作两性离子是不同物种，必须选择主分支；
- `F3_021`：确认二盐酸盐化学计量、质子化位置与 pD；
- 任一不确定状态都可建立 sensitivity branch，但不能混成一个平均谱。

## F3_023：身份冲突

当前 registry 采用 `3-(2-hydroxybenzamido)-L-tyrosine zwitterion` 结构。旧资料存在命名/结构冲突。必须由样品所有者根据原始合成记录、批次和分析证据确认；确认前所有计算保持 `IDENTITY_UNCONFIRMED`，不能用于最终结构结论。

## 一个显式水的含义

`+2 DMSO +1 H2O` 只是在给定局部氢键 motif 下的敏感性模型。它不表示样品中只有一个水，也不提供真实水含量、活度、聚集数或交换动力学。不同固定组成并排报告；只有获得独立的平衡/人口模型后才考虑跨组成组合。

## LVK 附录不进入 NMR 门控

先前报告中的 embedding、波形、时间和天空证据属于另一项 LVK 候选对分析。它与分子 NMR 的势能面、屏蔽常数、J 耦合和溶液构象没有科学联系。本包不复制其评分，也不允许它影响 NMR 候选、权重或方法选择。
