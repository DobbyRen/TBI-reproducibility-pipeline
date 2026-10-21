# 文章 I · 统一预处理与七级评估的 R 脚本骨架

对应手稿 **Methods B 部分 §2.5–§2.10**，以及预注册方案中的冻结参数。
目标是：所有"可能事后被质疑改过"的判据都写在代码常量里，跑出来的每一张表都能追溯到那一行。

---

## 1. 文件与设计文档的对应

| 手稿章节 | 文件 | 主要函数 |
|---|---|---|
| §2.5 统一预处理 | `03_preprocess.R` | `audit_scale()` `collapse_probes()` `log2_if_needed()` `filter_by_expression()` `maybe_combat()` `preprocess_unified()` |
| §2.6 统一重算 | `04_de.R` | `run_de_array()` `run_de_rnaseq()` `run_de_voom()` `run_de()` `threshold_grid()` `rebuild_modules()` |
| §2.7 内部稳定性 | `05_stability.R` | `patient_level_split()` `dual_fold_de()` `size_matched_draw()` `negctrl_compare()` |
| §2.8 七级度量 L1–L6 | `06_seven_levels.R` | `level1_burden()` `level2_membership()` `level3_robustness()` `level3_topk()` `level4_concordance()` `level5_stratified()` `level6_pathway()` |
| §2.9 判别性能与操作点（L7） | `07_model_transfer.R` | `run_l7()` `lock_operating_point()` `evaluate_locked()` `summarise_l7()` |
| §2.10 影响因素 | `08_determinants.R` | `fit_determinant()` `build_predictors()` `run_determinants()` |
| 冻结参数 | `00_config.R` | `FROZEN` `log_frozen()` `seed_for()` `assert_primary_threshold()` |
| 工具 / 取数 / 编排 | `01_utils.R` `02_io.R` `run_all.R` | — |

**关键设计：`07_model_transfer.R` 里 `lock_operating_point()` 的签名只接受训练数据。**
这不是风格问题，而是把「在接触外部数据之前锁定阈值」这条预注册承诺**写进了函数接口**——物理上拿不到外部数据，也就不可能事后回调阈值。

---

## 2. 环境安装

本机没有 R（已确认 `Rscript` / `Rterm` 均不存在）。先装 R ≥ 4.3（建议 4.4），再在 R 里执行：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# 核心（必须）
BiocManager::install(c("GEOquery", "limma", "DESeq2", "sva", "Biobase", "AnnotationDbi"))
install.packages(c("glmnet"))

# 可选（缺失时脚本会跳过或降级，不会报错中断）
BiocManager::install(c("WGCNA", "clusterProfiler", "biomaRt"))
install.packages(c("ranger", "xgboost", "betareg", "sandwich", "lmtest"))
```

降级路径已写好：DESeq2 缺失 → 自动退化为 `limma-voom`；ranger/xgboost 缺失 → 跳过敏感性模型；sva 缺失 → 跳过 ComBat 并在 audit 表中记录原因。

---

## 3. 目录结构

```
TBI_repro_pipeline/
├── 00_config.R            冻结参数（预注册后禁止修改）
├── 01_utils.R             集合度量、自助法、AUC/PR-AUC/校准
├── 02_io.R                清单读取、GEO 下载缓存、合成数据
├── 03_preprocess.R        统一预处理五步
├── 04_de.R                统一差异表达
├── 05_stability.R         患者水平双层拆分 + 负对照
├── 06_seven_levels.R      L1–L6
├── 07_model_transfer.R    L7 与锁定操作点
├── 08_determinants.R      影响因素（beta 回归 / 两部分模型）
├── run_all.R              端到端编排
├── check_syntax.py        语法自检（无需 R）
├── data/
│   └── datasets_template.csv
└── outputs/               运行后生成
```

---

## 4. 三种运行方式

```r
setwd("TBI_repro_pipeline 所在目录")

# (a) 合成数据自检：不联网、不依赖任何 Bioconductor 包也能跑通大部分环节
source("run_all.R", encoding = "UTF-8")     # 默认 synthetic

# (b) 真实分析
#     先核对 data/datasets_template.csv，再运行
#     Rscript run_all.R geo
```

**建议第一次一定先跑 (a)。** 它会在 `outputs/` 下生成全部表格，证明链路通畅，你再根据报错补装缺失包。

语法自检（不需要 R）：

```bash
python check_syntax.py
```

---

## 5. 与预注册参数的对应（哪些被冻结）

`00_config.R` 中的 `FROZEN` 列表逐项对应预注册 §5 与 §7：

- 主阈值 `FDR < 0.05 且 |log2FC| ≥ 0.5`
- 敏感性网格 `FDR ∈ {0.01, 0.05, 0.10}`
- 过滤规则：`≥20% 样本高于固定下限`
- 探针聚合：主分析 `maxmean`，敏感性 `mean`
- 批次：**仅在数据集显式提供批次变量时**才做 ComBat，且保护分组变量
- 拆分：患者水平、`prop = 0.5`、分层、固定种子
- 操作点：主分析 Youden，敏感性固定敏感度 90%
- 失效判据：锁定操作点下 sens 或 spec < 0.10

`log_frozen()` 会把 `00_config.R` 的 **MD5 与全部参数**写进 `outputs/frozen_params.log`。
投稿时把这个日志作为补充材料，就能证明「跑结果时用的就是这一套参数」。

**版本关系（2026-10-19 补记）：** OSF 归档的 `repo_skeleton.zip` 是注册时点的**冻结快照**，不可改动；本目录（及将来的 GitHub 仓库）是活代码库。两者分工——归档快照证明「注册时点的代码」，活仓库 + `frozen_params.log` 证明「产出结果的代码」。只要 `00_config.R` 不改，两者的分析口径永远一致。

---

## 6. 需要你人工完成的部分

| 项目 | 说明 |
|---|---|
| `data/datasets_template.csv` | 6 个数据集的分组字段名、对照/损伤水平取值**必须逐个核对**。我填的是占位值，`standardise_pheno()` 用模糊匹配定位水平，匹配失败会直接报错并列出实际取值，照着改即可 |
| `claimed_genes` 列 | 从提取表填入原始研究点名的 hub 基因，分号分隔。这是 L5 的分子 |
| `reported_de_n` 列 | 原始研究**自己报告**的 DEG 数，L1 的比值要用 |
| `data/gmt/*.gmt` | L6 需要基因集。从 MSigDB 下载 Hallmark / KEGG / GO_BP 的 `.gmt` 放进去，否则 L6 会跳过 |
| 单细胞数据集 | GSE269748 / GSE160763 需**先做 pseudobulk 聚合**（按细胞类型或样本），不能直接进本流水线 |
| 方向校准 | 见下 |

---

## 7. 三处必须你拍板的点

**① 阈值网格与预注册不完全一致 —— 已定处置口径，照此执行。**
注册协议 PDF 写的是 `|log2FC| ∈ {0.25, 1.00}`（敏感性），而**随注册归档的代码**（`repo_skeleton.zip` 内 `00_config.R`）实现的是并集 `{0, 0.25, 0.50, 1.00}`。协议 PDF 与归档代码都已被注册锁定、均不可改动；处置口径为**以归档代码为准**（更具体的操作文档优先），主分析阈值 FDR<0.05 且 |log2FC|≥0.5 三处完全一致，多出的 0 点是无效应参考点。手稿 §3.2 L3 表已改为与代码逐字一致，该差异已写进手稿 Results 的 *Deviations from the registered plan* 小节申报。**不要在提交后偷偷把网格改成其中任何一份。**

**② L5 的方向判据目前是硬编码的 `sign(logFC) == +1`。**
意思是"在外部数据中同样上调才算确认"。但各数据集的对比方向（损伤 vs 对照，还是对照 vs 损伤）未必统一。
正式跑之前必须按各数据集的实际对比方向统一校准，否则方向判据会把真信号判成未确认。`run_all.R` 结尾也会打印这条提醒。

**③ 负对照的定义方式要在论文里写清楚。**
`size_matched_draw()` 从**全队列排名**中随机抽同样数量的基因——这是 REDE 的原做法，目的是区分"拆分策略的贡献"和"集合变小的贡献"。审稿人一定会问这个，所以 §2.7 必须把"随机抽样的来源集合"写明确。

---

## 8. 与文章 II 的接口

`outputs/internal_stability.csv` 里的 `n_concordant` 基因集，就是文章 II 要用的**可复现特征白名单**。
`05_stability.R` 中 `dual_fold_de()` 返回的 `internal_set` 已按"两折均达主阈值且方向一致"筛选，可直接作为文章 II 的特征准入条件——这正是方案 C 里两篇文章的衔接点。

---

## 9. 已知边界

- 本仓库**未在本机执行过**（本机无 R）。`check_syntax.py` 只做括号/引号配平检查，不等于 R 能跑通。第一次运行请务必用 synthetic 模式。
- `06_seven_levels.R` 的 L6 用的是自实现的超几何 ORA；论文正式版建议改用 `clusterProfiler` 并保留自实现版本作为交叉校验（两者数学等价）。
- `08_determinants.R` 在合成数据模式下不运行（需要真实的研究层面特征）。
