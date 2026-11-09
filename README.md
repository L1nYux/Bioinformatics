# 已知母蛋白中的片段定位：三种比对策略的比较

从 10 条人类蛋白质中截取已知位置的片段，改变片段长度与氨基酸替换比例，比较严格全局、局部和半全局比对找回原始起止位置的能力。论文、输入序列、逐条预测和绘图程序均在本仓库中。[结课论文全文](结课论文.md)

## 实验结果

每条母蛋白抽取 3 个起点，每处截取 20、40、80 个残基，各施加 0%、10%、30% 三档替换。三种方法使用相同的 270 个查询，以各自返回的第一个最优解计分。

| 方法 | 两端均正确 | 精确定位率 | 平均覆盖率 | 平均端点误差（残基） |
|---|---:|---:|---:|---:|
| 严格全局 | 206/270 | 76.30% | 100.00% | 12.467 |
| 局部 | 215/270 | 79.63% | 99.13% | 0.154 |
| 半全局 | 268/270 | 99.26% | 100.00% | 0.006 |

以乳酸脱氢酶片段 `P00338_012_20_30` 为例，真实区间是 `[12,32)`。严格全局比对覆盖了全部 20 个残基，却将起点放在 0；局部比对从 13 开始；半全局比对找回了两个端点。半全局的两次首解偏差都存在同分的正确区间。结论限于已知母蛋白、只含人工替换的本组样本。[逐条预测](results/results.csv) · [等分区间](results/ambiguous.csv)

## 复现实验

在 Windows PowerShell 中，从解压后的项目根目录执行。以下以安装了 Python 启动器的机器为例；若使用独立 Python 的绝对路径，可用该解释器替代第一行的 `py -3.13`。

```powershell
py -3.13 -m venv .venv
& .\.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
& .\.venv\Scripts\python.exe fragment.py data/raw/P69905.fasta examples/query.fasta --mode semiglobal
& .\.venv\Scripts\python.exe -m pytest -q tests
& .\.venv\Scripts\python.exe experiment.py
& .\.venv\Scripts\python.exe analyze.py
```

示例输出区间为 `[20,40)`，报告坐标为第 21–40 位；查询 20 个残基均参与配对。单次模式也接受 `global`、`local`。输入为单条 FASTA，仅接受 20 种标准氨基酸。覆盖率为双方均有残基的比对列数除以查询长度。

## 复查结果

| 文件 | 用途 |
|---|---|
| [study.json](study.json) | 锁定评分、长度、替换比例、种子和起点数 |
| [study_manifest.csv](data/study_manifest.csv) | 10 条原始 FASTA 的来源、版本、长度与 SHA-256 |
| [queries.csv](results/queries.csv) | 查询序列、母蛋白、真值及从零计数的替换位置 |
| [results.csv](results/results.csv) | 每个查询每种方法一行，首个最优解的位置、覆盖率、得分和误差 |
| [per_protein.csv](results/per_protein.csv) / [summary.csv](results/summary.csv) | 各母蛋白的指标及跨蛋白汇总 |
| [ambiguous.csv](results/ambiguous.csv) | 存在多个等分端点区间的查询和模式 |

实验在 Python 3.13.5 / Windows 上完成。若要保留既有结果，给 `experiment.py` 和 `analyze.py` 同时传入 `--output 复现结果`。绘图读取指定目录中的 CSV。在固定版本的 Biopython 下取第一个最优解；复核阶段枚举等分解，本数据集已完成枚举。源码包附 `SHA256SUMS.txt`，哈希基于原始文件字节。每次正式实验读取本地快照。

## 署名和许可

本项目程序代码采用 [MIT 许可](LICENSE)。母蛋白数据来自 UniProt Consortium 的 UniProtKB/Swiss-Prot，按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 提供。原始条目链接逐条列于数据清单。项目对这些序列截取片段、人工替换部分氨基酸并计算预测；修改记录保存在查询表中。[数据方说明](https://www.expasy.org/resources/uniprotkb-swiss-prot)
