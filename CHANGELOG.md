# CHANGELOG —— 注册后代码变更台账

**对应注册：** OSF `10.17605/OSF.IO/GFDJR`（2026-10-18 提交）
**归档快照：** `osf_attachments/repo_skeleton.zip`（注册时点冻结，不可改动）
**本目录：** 活代码库。此处逐条记录活库相对归档快照的改动，投稿时据此撰写 *Deviations from the registered plan*。

**红线：** `00_config.R` 一律不改。当前 MD5 = `8562307b993ac2470f764db065994a00`，与归档快照内同名文件逐字节一致（2026-10-21 已校验）。

---

## 2026-10-21

### 1. 修复：`run_all.R` 无法从命令行启动（bug）
- **现象：** `Rscript run_all.R synthetic` 立即报 `Error in ensure_dirs(...) : could not find function "ensure_dirs"`。
- **原因：** 第 22 行在 `source("01_utils.R")` **之前**就调用了 `01_utils.R` 里定义的 `ensure_dirs()`。
- **改动：** 把建目录的调用挪到 source 循环之后（仅调换顺序，未改任何逻辑）。
- **影响：** 无。不触及 `00_config.R`，不改变任何参数、任何统计量、任何输出表的定义。此修复只是让注册时写好的代码真正能跑起来。

### 2. 新增：两个环境辅助脚本（不参与分析）
- `install_deps.R` —— 幂等地安装 12 个 Bioconductor 包 + 9 个 CRAN 包，末尾打印核心/可选包版本清单。
- `make_gmt.R` —— 用 `msigdbr` 生成 Hallmark / KEGG / GO_BP 三个 `.gmt` 到 `data/gmt/`（免 MSigDB 登录）。
- 二者都是**环境准备脚本**，不进入数据分析链路，不产生任何论文结果。

### 3. 环境事实（写进论文 Methods / 补充材料）
- R 4.6.1 (2026-06-24)，平台 `x86_64-apple-darwin20`；Bioconductor 3.23。
- 包库路径：`/Users/renqingqing/R/x86_64/4.6/library`（macOS 默认用户库目录不存在时 R 会静默降级到系统库，已在 `~/.Renviron` 固定）。
- 每次运行自动生成 `outputs/frozen_params.log`（含 `00_config.R` 的 MD5 与全部冻结参数）与 `outputs/session_info.txt`。

### 4. 首次自检状态
- `synthetic` 模式跑通，产出 L1–L5 全部表格。
- **L6 未执行**：`data/gmt/` 为空（基因集包下载中）。
- **L7 未执行**：合成数据双层拆分后交集基因 < 3，代码按设计跳过。→ 待补强信号冒烟测试。

### 5. 待决定（尚未改动代码）
跨物种基因标识：6 个数据集跨 human / rat / mouse 三物种，而 `07_model_transfer.R` 里的 `ortholog_map()`（biomaRt 一对一同源映射）**从未被调用**。是否接入、如何接入，待定后另记一条。
