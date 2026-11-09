# 蛋白质片段定位实验

这是一套用于生物信息学课程的小型 Python 工具。输入一条已知母蛋白和一个查询片段，可以查看三种序列比对方法给出的位置；使用随包数据，可以从头重建论文中的 270 个查询、810 条预测、分组表和两幅图。第一份可分享源码包版本为 0.1.0。

## 从哪里开始

先读[结课论文](结课论文.md)第 3 节的 P00338 案例，里面同一条查询有三个不同的预测起点。想理解成因，从单次命令开始；想核对论文数字，运行批量命令。安装时需要下载依赖，之后的比对、测试和绘图可以离线完成。

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

## 自己练习 AI 辅助编程

1. 先独立预测：给母蛋白两端各加一个残基，哪种模式的分数应保持不变？向 AI 提出“构造只延长母蛋白、再交换长短序列的反例”，然后对照 [test_fragment.py](tests/test_fragment.py) 的断言。
2. 在 [queries.csv](results/queries.csv) 中定位 `P00338_012_20_30`，先手算三种方法的端点误差，再把自己的计算和结果行交给 AI 核对。你应该得到 6、0.5 和 0 个残基。
3. 在 [ambiguous.csv](results/ambiguous.csv) 中找一个首解失败但存在精确解的例子。给 AI 的任务是“解释两个最优区间为什么都可以成立，以及软件应该向只有查询序列的用户显示什么”，再用实际区间审阅它的回答。

## 署名和许可

本项目程序代码采用 [MIT 许可](LICENSE)。母蛋白数据来自 UniProt Consortium 的 UniProtKB/Swiss-Prot，按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 提供。原始条目链接逐条列于数据清单。项目对这些序列截取片段、人工替换部分氨基酸并计算预测；修改记录保存在查询表中。[数据方说明](https://www.expasy.org/resources/uniprotkb-swiss-prot)
