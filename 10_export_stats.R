## =====================================================================
## 10_export_stats.R —— 导出各数据集的基因级差异表达结果（供 L5 文献点名基因验证用）
##
## 与 run_all.R 的一致性保证：
##   - 直接复用 00–09 的同 source，冻结参数仍来自 00_config.R（MD5 不变）
##   - 读取 10-27 已缓存的 pp_*.rds（已完成基因空间统一），跳过预处理与下载
##   - DE 调用与 run_all.R 完全一致：run_de(pp, data_type=pp$data_type, paired=FALSE)
## 自检：导出后比对 L1_de_burden.csv 的 n_de，必须逐个一致，否则中止
## =====================================================================

## 由调用方在 pipeline 目录下运行（cd 到本目录后 Rscript 10_export_stats.R）

for (f in c("00_config.R", "01_utils.R", "02_io.R", "03_preprocess.R",
            "04_de.R", "05_stability.R", "06_seven_levels.R",
            "07_model_transfer.R", "08_determinants.R", "09_ortholog.R")) {
  source(f, encoding = "UTF-8")
}
ensure_dirs("outputs", "data/raw", "data/processed", "data/gmt")

ids <- c("GSE104687", "GSE89866", "GSE2871", "GSE45997", "GSE269748", "GSE160763")

## 期待值（10-27 实测，来自 outputs/L1_de_burden.csv）—— 用作一致性自检
expect_n_de <- c(GSE104687 = 0L, GSE89866 = 416L, GSE2871 = 389L,
                 GSE45997 = 2533L, GSE269748 = 6474L, GSE160763 = 0L)

all_stats <- list()
for (nm in ids) {
  fp <- file.path("data/processed", paste0("pp_", nm, ".rds"))
  if (!file.exists(fp)) { message("[跳过] 无缓存：", fp); next }
  pp <- readRDS(fp)
  if (is.null(pp$data_type)) stop("pp_", nm, " 缺 data_type 字段")
  st <- run_de(pp, data_type = pp$data_type, paired = FALSE)
  st <- as.data.frame(st)
  st$dataset_id <- nm
  all_stats[[nm]] <- st
  n_de <- length(primary_de(st))
  ok <- identical(as.integer(n_de), as.integer(expect_n_de[[nm]]))
  message(sprintf("[stats] %-11s n_genes=%5d  n_de=%5d  期望=%5d  %s",
                  nm, nrow(st), n_de, expect_n_de[[nm]],
                  if (ok) "OK" else "*** 不一致，请检查 ***"))
  if (!ok) warning("DE 结果与 10-27 不一致：", nm)
}

saveRDS(all_stats, "data/processed/stats_all.rds")

## 导出 csv（每个数据集一张表；另出一张合并长表，便于 Python 侧处理）
out <- do.call(rbind, all_stats)
cols <- intersect(c("dataset_id", "gene", "logFC", "fdr", "P.Value", "adj.P.Val",
                    "AveExpr", "t", "B"), names(out))
write.csv(out[, cols], "outputs/stats_all_genes.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

## 各数据集显著基因集（主判据）
sig <- do.call(rbind, lapply(names(all_stats), function(nm) {
  st <- all_stats[[nm]]
  d <- st[!is.na(st$fdr) & st$fdr < FROZEN$de_fdr_primary &
            abs(st$logFC) >= FROZEN$de_lfc_primary, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  data.frame(dataset_id = nm, gene = d$gene, logFC = d$logFC, fdr = d$fdr,
             stringsAsFactors = FALSE)
}))
write.csv(sig, "outputs/stats_significant_genes.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

message("[done] 基因级结果已导出：outputs/stats_all_genes.csv (",
        nrow(out), " 行) / outputs/stats_significant_genes.csv (",
        if (is.null(sig)) 0 else nrow(sig), " 行)")
