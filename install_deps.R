# ---------------------------------------------------------------
# install_deps.R —— 文章 I 流水线依赖安装（一次性）
# 用法：Rscript install_deps.R
# 幂等：已装的包会跳过，可反复运行；失败时直接重跑即可
# ---------------------------------------------------------------

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

CRAN_PKGS <- c(
  "glmnet",                                    # L7 主模型（硬依赖）
  "ranger", "randomForest", "xgboost", "gbm",  # L7 敏感性模型（可缺）
  "betareg", "lmtest", "sandwich",             # 08_determinants.R
  "WGCNA"                                      # 模块重建（原研究用过才触发）
)

BIOC_PKGS <- c(
  "GEOquery", "limma", "DESeq2", "sva", "Biobase", "AnnotationDbi",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "biomaRt",
  "impute", "preprocessCore"   # WGCNA 的依赖在 Bioconductor，先从这里装，否则 WGCNA 会失败
)

have <- function(p) suppressWarnings(requireNamespace(p, quietly = TRUE))

cat("[install] R:", R.version.string, "\n")
cat("[install] Bioconductor:", as.character(BiocManager::version()), "\n\n")

for (p in BIOC_PKGS) {
  if (have(p)) { cat("[skip ]", p, "already installed\n"); next }
  cat("[bIOC ] installing", p, "...\n")
  ok <- try(BiocManager::install(p, ask = FALSE, update = FALSE), silent = TRUE)
  cat(if (have(p)) "[ ok  ] " else "[FAIL ] ", p, "\n", sep = "")
}

for (p in CRAN_PKGS) {
  if (have(p)) { cat("[skip ]", p, "already installed\n"); next }
  cat("[cran ] installing", p, "...\n")
  try(install.packages(p, quiet = TRUE), silent = TRUE)
  cat(if (have(p)) "[ ok  ] " else "[FAIL ] ", p, "\n", sep = "")
}

# ---- 验证：缺哪个一眼看见 ----
CORE <- c("GEOquery", "limma", "DESeq2", "sva", "Biobase",
          "AnnotationDbi", "biomaRt", "glmnet")
OPT  <- c("betareg", "lmtest", "sandwich", "WGCNA",
          "ranger", "randomForest", "xgboost", "gbm")

cat("\n[verify] --- 核心包 ---\n")
for (p in CORE) cat(sprintf("  %-14s %s\n", p,
  if (have(p)) as.character(packageVersion(p)) else "*** MISSING ***"))
cat("[verify] --- 可选包 ---\n")
for (p in OPT) cat(sprintf("  %-14s %s\n", p,
  if (have(p)) as.character(packageVersion(p)) else "(missing, 可降级)"))

miss <- CORE[!vapply(CORE, have, logical(1))]
if (length(miss)) {
  cat("\n[verify] 核心包缺失，流水线会停：", paste(miss, collapse = ", "), "\n")
  quit(status = 1)
} else {
  cat("\n[verify] 全部核心包就位。下一步：Rscript run_all.R synthetic\n")
}
