# CHANGELOG —— 注册后代码变更台账

**对应注册：** OSF `10.17605/OSF.IO/GFDJR`（2026-10-18 提交）
**归档快照：** `osf_attachments/repo_skeleton.zip`（注册时点冻结，不可改动）
**本目录：** 活代码库。此处逐条记录活库相对归档快照的改动，投稿时据此撰写 *Deviations from the registered plan*。

**红线：** `00_config.R` 一律不改。当前 MD5 = `8562307b993ac2470f764db065994a00`，与归档快照内同名文件逐字节一致（2026-10-21 提交前再次校验通过）。

---

## 2026-10-21

### 1. 修复：`run_all.R` 无法从命令行启动（bug）
- **现象：** `Rscript run_all.R synthetic` 立即报 `could not find function "ensure_dirs"`。
- **原因：** 第 22 行在 `source("01_utils.R")` **之前**就调用了 `01_utils.R` 里定义的 `ensure_dirs()`。
- **改动：** 把建目录的调用挪到 source 循环之后（仅调换顺序，未改任何逻辑）。
- **影响：** 无。不触及 `00_config.R`，不改变任何参数、统计量或输出表定义；只是让注册时写好的代码真正能跑起来。

### 2. 修复：L6 在 fdr 含 NA 时崩溃（bug）
- **现象：** 基因集就位后，L6 报 `Error in if (sum(eA$sig) > 0) ... argument is of length zero`。
- **原因：** `eA$sig <- eA$fdr < FROZEN$enrich_fdr` 在 `fdr` 为 NA 时得到 NA，`sum()` 随之返回 NA，`if (NA > 0)` 直接报错。**这条路径在真实数据上必然会触发**（某通路在两数据集公共基因空间内为空时，超几何 p 值为 NA）。
- **改动：** `sig` 显式判 `!is.na(fdr) & fdr < 阈值`；`replicated` 增加 `direction` 的 NA 保护；求和一律 `na.rm = TRUE`。
- **影响：** 只影响 NA 的处理口径（NA → 判为非显著），不改变任何有效数值的算式。

### 3. 新增：`09_ortholog.R` —— 跨物种基因标识统一（方案 A，**本次定案**）
- **背景：** 6 个数据集跨 human / rat / mouse。同一基因的 symbol 大小写不同（人 `TNF` / 鼠 `Tnf`），直接按字符串求交集会把跨物种对的重叠人为压到 0，污染 L4/L5/L6/L7。
- **方案 A 定义：** 以**人**为目标基因空间；非人数据集经 Ensembl **一对一同源**关系重命名；只保留严格一对一（两个方向都不允许重复），其余丢弃并记入审计表。
- **挂载位置：** 在**差异表达之前**、作用于表达矩阵行名（不能只改 DE 列表——L5/L7 还要用矩阵，两侧必须落在同一基因空间）。
- **失败策略（重要）：** 统一失败的数据集**从全部跨数据集层级中剔除**，绝不静默退回原基因空间，避免把被污染的数字当结果。
- **实现为何绕开 `biomaRt::getLDS`（实测，勿改回）：**
  1. `getLDS` 双 mart 查询返回 **HTTP 405 Method Not Allowed**；
  2. biomaRt 会把请求重定向到它 pin 住的归档站（本次 `jun2026.archive.ensembl.org`），且 curl 超时仅约 10 s，必超时；`biomaRt 2.68.0` 的 `useMart()` 无 `ensemblRedirect` 参数可关。
  → 改为直接请求 BioMart `martservice` 端点（跟随 308 重定向、超时 240 s），一次成功（16–25 s）。
- **产出与可复现性：** 映射表缓存到 `data/processed/ortholog_<from>_to_<to>.rds`（入库），并导出 `outputs/ortholog_map_*.csv`；审计表 `outputs/ortholog_audit.csv` 记录每数据集的物种、状态、丢弃数、折叠数与 Ensembl release。缓存命中后不再联网，结论不受 Ensembl 版本漂移影响。
- **实测规模：** mouse→human 严格一对一 **17,041** 对；rat→human **16,444** 对；Ensembl release **jun2026**。

### 4. 新增：`test_smoke_crossspecies.R` —— 强信号冒烟测试
覆盖两条此前**从未被执行过**的路径：
- **跨物种统一**：用真实同源表构造一对小鼠/人合成数据集，共享同一批同源对与同一效应向量。结果：DE 集合重叠 **134 / 150**，Jaccard **0.812**（不做映射则必为 0）。
- **L7 锁定操作点迁移**：youden 规则下外部 AUC 1.000、sens 0.950 / spec 1.000，锁定点**未**失效；fixed_sens90 下 sens 恰好 0.900。
- 同时验证 L6（Hallmark 50 个集合）可计算。
- **注意：** 首次实现时曾出现 `sens = 0 / spec = 1`，排查后确认是**模拟本身的缺陷**（两侧效应向量独立抽取，方向可能相反），非代码问题；改为共用效应向量后恢复正常。此结论须如实写进测试脚本注释，避免日后误判为 L7 缺陷。

### 5. 修改：`make_gmt.R` 改用 Broad 开放镜像
- **原路径失败：** MSigDB 官网下载**需登录**，未登录时返回的是 4425 字节的 HTML 登录页（HTTP 200，极易误判为成功）；`msigdbr` 的 zenodo 源实测约 10 KB/s，下不动。
- **现路径：** `data.broadinstitute.org/.../msigdb.v7.5.1.symbols.gmt`（28 MB，实测 5.7 MB/s，无需登录，一次拿到全部集合）。
- **产出：** `data/gmt/msigdb_{Hallmark,KEGG,GOBP}_v7.5.1.gmt`（50 / 186 / 7658 个集合）+ `README_provenance.md`（记录 URL、md5、拆分时间）。
- **需申报的偏离：** 注册协议未指定 MSigDB 版本，此处用 **v7.5.1**（Broad 开放镜像可提供的版本），须在 Methods 写明版本与获取途径。

### 6. 环境事实（写进 Methods / 补充材料）
- R 4.6.1 (2026-06-24)，平台 `x86_64-apple-darwin20`；Bioconductor 3.23；biomaRt 2.68.0。
- 包库路径 `/Users/renqingqing/R/x86_64/4.6/library`（macOS 默认用户库不存在时 R 会静默降级到系统库，已在 `~/.Renviron` 固定 `R_LIBS_USER`）。
- 每次运行自动生成 `outputs/frozen_params.log`（含 `00_config.R` 的 MD5 与全部冻结参数）与 `outputs/session_info.txt`。

### 7. 已知方法学注记（须写进论文 Limitations）
- `cv.glmnet` 在"每折样本 < 10"时会把 CV 准则从 `auc` 自动降级为 `deviance`（本次运行触发 30 次该警告）。这意味着小样本数据集上，L7 的 λ 选择实际以 deviance 为准而非 AUC，需在 Limitations 中说明。
- L5 的方向判据目前固定为 `sign(logFC) == +1`，**正式出结果前必须按各数据集实际对比方向统一校准**（`run_all.R` 收尾处也有该提醒）。

### 8. 待办
- GitHub 仓库推送：本地已 commit（`16fd051`），远端待用户在 GitHub 建空仓库后 push（须用 Personal Access Token，本机无 SSH 密钥）。
- `data/gmt` 中 GO_BP 集合数达 7658，正式运行时 L6 的组合数较大，需关注运行时间。
