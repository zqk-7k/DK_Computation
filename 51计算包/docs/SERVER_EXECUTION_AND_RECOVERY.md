# 服务器执行顺序与稳健恢复

本文给出科学流程的执行约束。具体 Slurm/PBS 资源值以 `hpc/` 中站点 profile 为准；不能把一个集群的 CPU/内存语义直接复制到另一个集群。

## 1. 上传前

1. 运行包审计，确认 51 个 ID，不含 `F3_005`；
2. 验证所有 JSON/CSV、结构哈希、输入 bundle hash 和 checksum inventory；
3. 确认没有账号密码、SSH 私钥、令牌或实验敏感信息；
4. 将 condition lock 中未确认字段保持 `PROVISIONAL`；
5. 只上传完整目录，不单独复制输入而丢失 config、gates 和 provenance。

## 2. 运行环境预检

- 固定 ORCA 6.1.1 可执行文件绝对路径并记录二进制哈希；
- 使用与 ORCA 匹配的 OpenMPI；调用串行 ORCA driver，不对 driver 执行 `mpirun`；
- 设置 `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`，除非特定步骤另有验证；
- 检查 basis、aux basis、NMR/J basis 和 relativistic route；
- 创建节点本地 scratch，确认可写和空间足够；
- 运行每种资源 profile 的最小 smoke test。

## 3. 顺序和依赖

严格遵循：condition/species lock → GOAT wave → 结构门/去重 → r2SCAN opt/freq → PBE0 opt/freq → QRRHO/人口冻结 → PBE0 shielding/reference → revDSD gate/reference → J gate → spin simulation → clean output → optional appearance。

下游作业只能读取已签署的上游 manifest 和 hash receipt。不能因某些结构先完成，就越过尚未终态的同波次任务。

## 4. GOAT 的安全早停

备用候选保持 HELD 或不提交。达到数量和覆盖门后只停止释放备用任务；不取消已运行任务。取消 running job 既浪费计算，也会引入完成时间偏差。

GOAT 通常不能像普通几何优化那样从任意内部迭代完全恢复。节点失败时保留日志和已输出结构，但重新提交应从经验证的起点或工作流定义的恢复点开始，并标记为新 attempt。

## 5. 内存预算原则

ORCA `%maxcore` 是近似每进程内存，不是作业总内存。保守门为：

`nprocs × maxcore_MB ≤ 0.70–0.75 × scheduler_requested_memory_MB`

其余内存留给积分、MPI、文件缓冲和程序峰值。增加 CPU 不一定缩短 MP2 response；更多 MPI ranks 可能提高复制内存和通信成本。revDSD 应先用少量代表结构做 4/8 ranks 的 scaling + memory pilot，再固定 profile。

## 6. 失败分类和动作

| 失败类型 | 先检查 | 正确动作 |
|---|---|---|
| `Read-only file system` / I/O | 节点磁盘、Slurm spool、stage-out | 改用已验证的 node-local scratch；原子化复制并校验哈希 |
| OpenMPI/OpenFabrics | ORCA/OpenMPI ABI、网络 transport | 加载匹配 MPI；仅在站点批准后禁用异常 transport |
| ORCA MP2 response 内存中止 | `%maxcore`、ranks、实际 MaxRSS | 降低 ranks/提高内存；先做 pilot，不把它误判为 Slurm OOM |
| scheduler OOM | cgroup/Slurm 日志 | 提高请求或降低进程和 maxcore，保留失败输出 |
| SCF 不收敛 | 波函数、电荷、自旋、初猜 | 独立 SCF 修复路线；不能只增加内存 |
| 虚频 | 振型、网格、优化阈值 | 沿模式重优化或提高网格；不静默取绝对值 |
| stale lock | scheduler job ID、心跳、节点状态 | 只有确认作业不存在且心跳超时后清理；记录操作者和时间 |
| walltime | checkpoint 可用性、阶段边界 | 在安全阶段拆分或从已批准几何重启；不要假装完整续算 |

## 7. revDSD 失败的保守恢复

对于先前 16 MPI × `%maxcore 6400` 在约 114 GB 节点上的失败，不能认为 16×6.4 GB 只占 102.4 GB 就安全；程序峰值和 MPI 开销可能超过剩余空间。优先使用站点资源 profile 中的 8-rank 高内存 pilot，并满足 70–75% 预算门、足够 scratch 和更长 walltime。若仍失败，再比较 4-rank 或高内存节点，不在同一输入上无分类地重复提交。

ORCA MP2-response 的 `.gbw` 不应被假定为可完整 checkpoint。恢复通常从已批准 XYZ/GBW 初猜重新运行 property，并保存前一 attempt 日志。

## 8. 文件落盘

所有作业在 scratch 中运行，结束后写临时目标目录，计算 SHA-256，再原子重命名为最终目录。每个 attempt 保存 stdout、stderr、ORCA input、geometry、版本、主机、资源、MaxRSS、walltime、退出码和哈希。只有 receipt 完整的结果才能进入下游 manifest。

## 9. 不应通过增加 CPU 解决的问题

错误物种、缺失 motif、错误参考、跨组成平均、交换峰、线宽/噪声和相对论路线都不是“多用 CPU”可以修复的。扩大资源前先确认瓶颈是可并行电子结构工作，而不是科学模型或内存带宽。
