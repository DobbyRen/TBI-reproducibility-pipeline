#!/usr/bin/env Rscript
## =====================================================================
## smoke_gse104687.R —— 单数据集集成冒烟：验证 loader 缓存 → 预处理 → DE
##
## 只跑 GSE104687 一个数据集，不改动 run_all.R、不改动 manifest。
## 目的：确认 build_geo_cache.R 造的缓存能被 02/03/04 正确消费，
##       并把过滤后的基因数、DE 数等实测值交给《执行日志》留痕。
##
## 用法（在 TBI_repro_pipeline 目录下）：
##   Rscript smoke_gse104687.R [region] [variant]
## =====================================================================
suppressMessages({
  for (f in c("00_config.R", "01_utils.R", "02_io.R", "03_preprocess.R",
              "04_de.R", "05_stability.R", "06_seven_levels.R",
              "07_model_transfer.R", "08_determinants.R", "09_ortholog.R")) {
    source(f, encoding = "UTF-8")
  }
})
ensure_dirs("outputs", "data/raw", "data/processed", "data/gmt")

args <- commandArgs(trailingOnly = TRUE)
REGION  <- if (length(args) >= 1) args[1] else "HIP"
VARIANT <- if (length(args) >= 2) args[2] else "rin"

GSE <- "GSE104687"
cache <- file.path("data", "raw", paste0(GSE, ".rds"))
if (!file.exists(cache)) stop("没有缓存，先跑 build_geo_cache.R。期望：", cache)

## 缓存里记的版本是否就是这次要的
esets  <- readRDS(cache)
parts  <- eset_to_parts(esets[[1]])
cat("[smoke] 缓存矩阵：", nrow(parts$mat), "基因 ×", ncol(parts$mat), "样本\n")
cat("[smoke] 脑区：", paste(sort(unique(parts$pheno$region)), collapse = "/"), "\n")

ph <- standardise_pheno(parts$pheno, "tbi_group", "control", "tbi",
                        patient_col = "donor_name", batch_col = NULL)
cat("[smoke] 表型：tbi =", sum(ph$group == 1L), "，control =", sum(ph$group == 0L),
    "；患者(供者)数 =", length(unique(ph$patient)), "\n")
cat("[smoke] 批次列：",
    if (all(is.na(ph$batch))) "NA（按协议：未证明混杂 → 不做 ComBat）" else "有值",
    "\n")

## 预处理（annotation = NA → 不聚合探针；本数据集行名已是 symbol）
pp <- preprocess_unified(parts$mat, ph, ann = get_annotation(parts$annot),
                         dataset_id = GSE)
cat("\n[smoke] ---- 预处理审计 ----\n")
print(t(pp$steps))

## 差异表达：FPKM（非计数）走 array 分支 = limma lmFit + eBayes
st <- run_de(pp, data_type = "array", paired = FALSE)
sig <- st$adj.P < FROZEN$de_fdr_primary & abs(st$logFC) >= FROZEN$de_lfc_primary
cat("\n[smoke] ---- DE 结果（主阈值 FDR<0.05 且 |log2FC|>=0.5）----\n")
cat("[smoke] 参与检验基因数 =", nrow(st), "；显著 =", sum(sig, na.rm = TRUE), "\n")
cat("[smoke] logFC 范围 =", round(range(st$logFC, na.rm = TRUE), 3), "\n")
if (any(sig, na.rm = TRUE)) {
  top <- st[order(st$adj.P), ][1:min(10, sum(sig, na.rm = TRUE)), ]
  cat("[smoke] 前 10 个（按 adj.P）：\n")
  print(top[, intersect(c("gene", "logFC", "P.Value", "adj.P"), names(top))])
}
cat("\n[smoke] 完成。region =", REGION, "（缓存实际内容以 build_geo_cache.R 参数为准）\n")
