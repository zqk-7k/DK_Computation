工作2-计算资源包-（3）
2026-10-13；64个18O产物的显式溶剂NMR流程

一、这次改变
命名：文件夹与ZIP均统一为“工作2-计算资源包-（3）”。后续修订继续编号。
原v1、v2和其压缩包没有删除、覆盖或改名；它们仅作为历史版本。
活动方法：以仓库51计算包的方法设计为基准，排除revDSD。
请以config/protocol.json和本说明为准，不要运行history/或reference_51/里的旧流程。

二、能否直接用2DMSO+1H2O？
可以。无需先完成0显式溶剂的计算，也无需先完成2DMSO/0水的计算。
工作模型是“产物+2个显式DMSO-d6+1个普通H2O”，同时放在DMSO连续介质中。
GOAT阶段用ALPB(DMSO)，DFT阶段用CPCM(DMSO)。
显式分子描述局部相互作用，连续介质描述外部总体溶剂环境，两者并存。
不可把它理解成先做连续介质、随后关闭连续介质改做显式溶剂。
三种组成是独立情景，没有串行依赖，不能跨情景归一化或混合玻尔兹曼权重。

保留的情景：
  DMSO_00_W00_continuum：0DMSO/0水；每模型3个起点；对照。
  DMSO_02_W00_primary：2DMSO/0水；每模型4个起点；无显式水模型。
  DMSO_02_W01_trace_water_sensitivity：2DMSO/1水；每模型7个起点；当前命令示例默认。
64产物对应66立体结构分支，另加TMS共67模型；201个组成分支、938个首批GOAT输入。
其中2DMSO+1水分支有469个首批输入（66产物模型+1TMS，各7个）。
2个DMSO不是测得的配位数或完整溶剂壳；1个水不是样品的定量含水浓度。
真实样品接近无水时，应保留不加水分支作为敏感性检查。不要声称含水模型必然更准确。

三、活动计算链（没有revDSD）
完整显式簇的GOAT/外部GFN2-xTB/ALPB采样
 -> 实际输出、原子顺序、立体/组成、氢键模式与采样审查
 -> r2SCAN-3c/CPCM预筛优化和频率
 -> PBE0-D4/def2-TZVPP/CPCM优化和频率
 -> 同一组成的qRRHO权重
 -> PBE0/pcSseg-3核磁屏蔽
 -> 同一情景、方法、版本和温度的TMS完整参比系综
 -> 两种独立J候选：PBE/pcJ-3与PBE0/pcJ-2
 -> 独立实验基准比较J，再进行方法选择
 -> 只汇总产物的1H与13C预测谱。

不将两种J方法的结果混合。含Br目标的J使用明确标记的Br/def2-TZVPP混合基组，
不能把普通轻元素基准直接当成含Br方法已验证。方法指纹对此作了区分。
采用所有通过检查且去重后的预筛构象进入PBE0，不做未经验证的激进截断。
与51包“至少保留95%且保留重要模式”的设计相比，这一保守选择可能更耗计算。

四、已实现与未完成必须区分
已实现：实际GOAT结果导入、显式簇DFT衔接、性质输入、J基准程序、目标/参比各自加权、
原子掩码、谱图汇总、收据与来源校验、失败停止、单任务/数组任务入口。
参数中保留18O精确质量和DMSO-d6的D质量；GOAT四列XYZ不写质量，后续内联坐标恢复。
不在谱图中加入显式溶剂的峰，也不把环境原子计入产物J耦合。

未完成：64个产物的实际GOAT/DFT、真实实验J基准、完整目标化合物精度验证。
本轮仅额外执行两个6原子水二聚体的真实ORCA输入测试：屏蔽和PBE/pcJ-3。
它们不验证DMSO溶剂壳、产物构象权重或实验准确度。日志在validation/real_smoke/。
不能把单元测试的合成夹具当成真实化合物结果。

仍需审查：F2_E11/F2_E22实际异构体、F4_D07中性体/盐型和溶剂、必要的互变异构态；
GOAT得到的簇是否保留合理组成，是否发生质子转移/解离，低频模式和采样是否充分。
当前298.15K是暂定值。更换真实温度等条件须建立新的协议版本，不混用结果。

五、入口与运行顺序
Python：本机用 py -3.13 -B -X utf8；Linux一般用python3。
需要RDKit、numpy、scipy、matplotlib及ORCA6.1.1；GOAT还需要兼容的外部xTB。
本次测试的Python依赖版本见requirements-tested.txt；服务器仍需单独检查ORCA/MPI/xTB运行环境。
本机缺少可用otool_xtb及其配套DLL；没有擅自安装，也没有启动64产物作业。

1. 看模型/状态：
   python workflow.py status

2. 仅生成2DMSO+1水的GOAT任务清单，不执行：
   python workflow.py goat-manifest --scenario DMSO_02_W01_trace_water_sensitivity

3. GOAT单作业执行入口见scripts/goat_portable.py --help。
   为避开外部程序对中文路径的兼容性差异，保留中文资源包名，但用ASCII工作目录。
   例如 --scratch-root E:/ORCA_work2_scratch。运行前先加--dry-run检查。
   服务器请改为自己的ASCII scratch路径，不要直接复制Windows盘符。
   只有显式run命令且确认模型假设后才执行；不会自动提交938个作业。

4. GOAT完成后生成审查草稿，再由操作者实际审查；草稿不是批准：
   python workflow.py draft-review --job-id JOB_ID --out reviews/review.json
   python workflow.py import-goat --job-id JOB_ID --review reviews/review.json
   在准备DFT前导入本轮全部已审查起点。不能计算后再悄悄增加样本并覆盖旧系综。

5. 为指定结构分支生成预筛输入（不执行）：
   python workflow.py prepare --species F2_E01_rep_R__DMSO
   未提供--scenario时，此版示例入口默认2DMSO+1水；TMS_DMSO也需相同情景单独准备。

6. 单个DFT作业运行：
   python workflow.py run --job JOB_RELATIVE_PATH --orca ORCA_EXECUTABLE --scratch-root ASCII_SCRATCH --acknowledge-models --dry-run
   确认后去掉--dry-run才运行。失败不会静默跳过；显式--retry保留旧日志。
   无收据的既有输出不能用--retry强行覆盖，须先单独审查；结果文件系统需支持同目录硬链接以原子提交收据。
   数组入口：先workflow.py manifest --stage STAGE，再run-task --manifest FILE --index INTEGER。
   每作业一个ORCA驱动；nprocs默认8，不能直接在外层再mpirun多个驱动。

7. 每一阶段全部检查通过后衔接：
   python workflow.py advance --species F2_E01_rep_R__DMSO --from-stage prescreen
   python workflow.py advance --species F2_E01_rep_R__DMSO --from-stage optfreq
   不会启动新作业，只生成下一步输入。遇到负频、原子改变、缺失输出即停止。

8. 真实结果完成后，可直接得到标记为“未验证预测”的表和谱：
   python workflow.py collect --species F2_E01_rep_R__DMSO --j-method PBE0_pcJ2 --exploratory
   这不会编造缺失数据。需要目标和匹配TMS全部真实结果；尚未计算时必须报错。
   输出为带来源的JSON/CSV、1H和13C谱SVG/PNG；谱图宽度是显示参数，不是预测弛豫。

9. 若要“经独立基准及收敛审查”的结果：
   提供真实实验J面板，workflow.py benchmark --panel PANEL.json --targets TARGET_IDS。
   两个候选方法均须有真实计算结果；报告不自动批准方法。
   正式collect还需--benchmark-gate与--convergence-gate，不能用程序伪造审查签名。
   协议说明与面板字段：docs/基准与收敛审查.txt。

六、自动化边界
机械步骤由程序完成；科学验收不预先盖PASS。
此版本采用“显式命令推进、失败即停”的阶段执行方式，不是无人值守的一键完成。
尚无自动扩充/放行后续GOAT波次的调度控制器，也不会自动取消运行中的任务。
采样扩展应使用独立版本/冻结波次记录；生产收敛必须有真实的来源绑定和稳定性证据。
任何计算预测都不能替代纯品身份、实验NMR和真实样品状态的确认。

七、历史与活动文件
活动配置：config/protocol.json；活动代码：workflow.py、clusterflow/、scripts/goat_portable.py。
初始构型：01_goat_explicit/；全部结构：structures_all_938_starts.sdf。
历史说明移至history/v2/；reference_51仅用于追溯，含被用户排除的方法不代表启用。
核心构型/来源文件保持原哈希。新版报告和压缩包校验由打包程序生成。
