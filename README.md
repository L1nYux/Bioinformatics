# 蛋白质片段定位：端部扰动与位置歧义

从 10 条人类蛋白质中截取已知位置的片段，比较严格全局、局部和半全局比对找回原始端点的能力。保留原 270 查询基线，新增五个种子、900 组匹配片段、3,600 个查询的端部干预实验，并实现有界的等分位置报告。[结课论文全文](结课论文.md) · [扩展实验报告](analysis/extension/report.md) · [项目检查与当前状态](docs/项目检查与改进_2026-11-10.md)

当前本地版本为 **0.2.0**。研究贡献集中于固定替换数量的配对设计、按蛋白分析和定位歧义的工程报告；经典比对方法采用 Biopython。结果限于已知母蛋白、均匀人工替换及当前评分，尚未评价跨物种、插入缺失或未知来源检索。

## 实验结果

每条母蛋白抽取 3 个起点，每处截取 20、40、80 个残基，各施加 0%、10%、30% 三档替换。三种方法使用相同的 270 个查询，以各自返回的第一个最优解计分。

| 方法 | 两端均正确 | 精确定位率 | 平均覆盖率 | 平均端点误差（残基） |
|---|---:|---:|---:|---:|
| 严格全局 | 206/270 | 76.30% | 100.00% | 12.467 |
| 局部 | 215/270 | 79.63% | 99.13% | 0.154 |
| 半全局 | 268/270 | 99.26% | 100.00% | 0.006 |

以乳酸脱氢酶片段 `P00338_012_20_30` 为例，真实区间是 `[12,32)`。严格全局比对覆盖了全部 20 个残基，却将起点放在 0；局部比对从 13 开始；半全局比对找回了两个端点。半全局的两次首解偏差都存在同分的正确区间。结论限于已知母蛋白、只含人工替换的本组样本。[逐条预测](results/results.csv) · [等分区间](results/ambiguous.csv)

扩展实验保持母片段、总替换数与共同的内部替换相同，将另外两个替换移至首端、末端或双端。每个位置组 900 个查询，首解结果如下。该对照单独保存，不能与原基线的采样分布直接混用。[配置](extension.json) · [逐条查询](analysis/extension/controlled_queries.csv)

| 替换位置 | 严格全局精确定位率 | 局部精确定位率 | 半全局精确定位率 |
|---|---:|---:|---:|
| 仅内部 | 84.78% | 96.00% | 99.67% |
| 首端 | 35.78% | 8.56% | 97.22% |
| 末端 | 24.22% | 8.67% | 98.56% |
| 双端 | 11.00% | 0.11% | 96.89% |

配对分析以母蛋白为簇计算差值和 bootstrap 描述区间，避免把相关片段当作独立蛋白。[完整统计与图表](analysis/extension/report.md)

## 复现实验

在 Windows PowerShell 中，从解压后的项目根目录执行。以下以安装了 Python 启动器的机器为例；若使用独立 Python 的绝对路径，可用该解释器替代第一行的 `py -3.13`。

```powershell
py -3.13 -m venv .venv
& .\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
& .\.venv\Scripts\python.exe fragment.py data/raw/P69905.fasta examples/query.fasta --mode semiglobal
& .\.venv\Scripts\python.exe -m pytest -q tests
& .\.venv\Scripts\python.exe experiment.py
& .\.venv\Scripts\python.exe analyze.py
& .\.venv\Scripts\python.exe extension.py
& .\.venv\Scripts\python.exe verify_project.py --reproduce
```

示例输出区间为 `[20,40)`，报告坐标为第 21–40 位；查询 20 个残基均参与配对。单次模式也接受 `global`、`local`。输入为单条 FASTA，仅接受 20 种标准氨基酸。覆盖率为双方均有残基的比对列数除以查询长度。

若只需检查当前交付，运行 `python verify_project.py`；`--reproduce` 会在项目内临时目录完整重跑两套实验，比较 CSV、图、扩展报告与运行清单的原始字节，并清理临时目录。复现需使用锁定依赖；依赖版本不符会直接报错。跨操作系统的图像字节一致性尚未验证。

## 查看唯一位置与等分候选

```powershell
& .\.venv\Scripts\python.exe fragment.py data/raw/P69905.fasta examples/query.fasta --mode semiglobal --json --max-alignments 1000
```

JSON 包含首解 `first`、去重端点 `intervals`、枚举状态及考察路径数。`unique` 表示完整枚举后的唯一端点，`ambiguous` 表示多个同分端点，`no_hit` 表示无配对，`incomplete` 表示路径上限内无法确认是否穷尽。即使上限内只看到一个位置，也会报告 `incomplete`，不会宣称唯一。默认文本输出也显示这些状态。

仅在 `unique` 时接受单一位置答复的规则不读取真实坐标。唯一最优端点仍可能错误：半全局在原基线接受 268/270 个查询且全部正确，在扩展双端组接受 877/900 个查询、其中 22 个错误。接受率、接受后错误率与片段覆盖率分别定义。离线结果列 `truth_in_observed_intervals` 仅供分析，不参与预测或接受决策。[基线统计](analysis/extension/baseline_summary.csv) · [扩展统计](analysis/extension/controlled_summary.csv)

## 复查结果

| 文件 | 用途 |
|---|---|
| [study.json](study.json) | 锁定评分、长度、替换比例、种子和起点数 |
| [study_manifest.csv](data/study_manifest.csv) | 10 条原始 FASTA 的来源、版本、长度与 SHA-256 |
| [queries.csv](results/queries.csv) | 查询序列、母蛋白、真值及从零计数的替换位置 |
| [results.csv](results/results.csv) | 每个查询每种方法一行，首个最优解的位置、覆盖率、得分和误差 |
| [per_protein.csv](results/per_protein.csv) / [summary.csv](results/summary.csv) | 各母蛋白的指标及跨蛋白汇总 |
| [ambiguous.csv](results/ambiguous.csv) | 存在多个等分端点区间的查询和模式 |
| [严格全局与局部分歧明细](analysis/global_local_discordance.csv) | 每个分歧查询一行，列出替换位置、真值、两种首解的端点与配对数，以及严格全局是否有端点正确的最优解 |
| [controlled_results.csv](analysis/extension/controlled_results.csv) | 10,800 次扩展比对的首解、去重候选、枚举状态与离线真值评价 |
| [paired_effects.csv](analysis/extension/paired_effects.csv) | 端部组相对内部组的配对成败、蛋白簇区间及逐种子范围 |
| [by_condition.csv](analysis/extension/by_condition.csv) / [by_seed.csv](analysis/extension/by_seed.csv) | 按长度、替换比例或随机种子分层 |
| [run_manifest.json](analysis/extension/run_manifest.json) | 扩展配置、程序哈希与实际依赖版本 |

实验在 Python 3.13.5 / Windows 上完成。若要保留既有结果，给 `experiment.py` 和 `analyze.py` 同时传入 `--output 复现结果`。绘图读取指定目录中的 CSV。在固定版本的 Biopython 下取第一个最优解；复核阶段枚举等分解，本数据集已完成枚举。源码包附 `SHA256SUMS.txt`，哈希基于原始文件字节。每次正式实验读取本地快照。

扩展复跑可使用 `extension.py --output 复现扩展`。该脚本仅从 `study.json` 读取评分参数，扩展种子与采样规则使用 `extension.json`。`release.py` 更新当前交付文件的 SHA-256 清单、运行检查，并生成 v0.2.0 源码与论文包。修改代码、配置或论文后，应先重跑受影响的实验，再执行 `release.py`。

分歧分析直接读取已有 CSV。下列命令按母蛋白和端部替换情况打印四种配对结果的计数，并另存 29 行分歧明细。查询的首残基位置为 0，末残基位置为长度减一。`local_only` 和 `global_only` 表示只有对应方法的首解精确定位；`global_truth_optimal` 表示至少一个严格全局最优解的端点与真值一致。基线实验及其原校验清单可按提交追溯。[基线版本 bf14ee1](https://github.com/L1nYux/Bioinformatics/tree/bf14ee1b433e5e3b43eab68d4aa8f1484734606a)

```powershell
& .\.venv\Scripts\python.exe compare_endpoints.py --input results --details 复现结果/global_local_discordance.csv
```

## 署名和许可

本项目程序代码采用 [MIT 许可](LICENSE)。母蛋白数据来自 UniProt Consortium 的 UniProtKB/Swiss-Prot，按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 提供。原始条目链接逐条列于数据清单。项目对这些序列截取片段、人工替换部分氨基酸并计算预测；修改记录保存在查询表中。[数据方说明](https://www.expasy.org/resources/uniprotkb-swiss-prot)
