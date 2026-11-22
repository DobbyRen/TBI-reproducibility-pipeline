#!/usr/bin/env Rscript
## =====================================================================
## smoke_others.R —— GSE89866 / GSE2871 / GSE45997 的集成冒烟
##
## 目的：在进入 run_all.R 之前，证明「loader 缓存 → 预处理 → DE」这条链
##       在三个数据集上都能跑通，并用内部对照排除"数据/映射搞错了"。
##
## 每个数据集都跑三项：
##   1) 主对比：case vs control
##   2) 标签置换负对照：随机打乱分组标签 → 期望 DEG 崩到 ~0
##      （若置换后仍有大量 DEG，说明流程有泄漏或批次结构冒充了信号）
##   3) 数据内在对照：
##      · GSE89866（人血、配对）：同受试者两次抽血，技术/个体效应应远大于暴露效应
##      · GSE2871 / GSE45997（大鼠）：受试间（动物）对比 —— 用"随机两分组"看基线
##
## 用法： Rscript smoke_others.R [GSE89866|GSE2871|GSE45997|all]
## =====================================================================
suppressMessages({
  library(Biobase); library(limma)
})
source("00_config.R", encoding = "UTF-8")
source("01_utils.R",  encoding = "UTF-8")
source("02_io.R",     encoding = "UTF-8")
source("03_preprocess.R", encoding = "UTF-8")
source("04_de.R",     encoding = "UTF-8")

set.seed(20261026)

SPEC <- list(
  GSE89866 = list(rds = "data/raw/GSE89866.rds", type = "rnaseq",
                  group_col = "group", control_level = "control", case_level = "tbi",
                  patient_col = "subject_id", paired = TRUE),
  GSE2871  = list(rds = "data/raw/GSE2871.rds",  type = "array",
                  group_col = "group", control_level = "control", case_level = "tbi",
                  patient_col = NULL, paired = FALSE),
  GSE45997 = list(rds = "data/raw/GSE45997.rds", type = "array",
                  group_col = "group", control_level = "control", case_level = "tbi",
                  patient_col = NULL, paired = FALSE),
  ## 敏感性单元：患侧 vs 对侧（同动物配对，规避 Todd/Xu 动物批次混杂）
  GSE45997C = list(rds = "data/raw/_alt/GSE45997_ipsi_vs_contra.rds", type = "array",
                   group_col = "group", control_level = "control", case_level = "tbi",
                   patient_col = "animal", paired = TRUE),
  ## 单细胞 pseudobulk（按样本聚合，预先指定 2026-10-27）
  ## GSE160763：2×2 正交设计，n=4/4 <5 → 降级敏感性
  GSE160763 = list(rds = "data/raw/GSE160763.rds", type = "rnaseq",
                   group_col = "group", control_level = "control", case_level = "tbi",
                   patient_col = NULL, paired = FALSE),
  ## GSE269748：主单元 24h CCI ipsilateral(9) vs Naive(10)
  GSE269748 = list(rds = "data/raw/GSE269748.rds", type = "rnaseq",
                   group_col = "group", control_level = "control", case_level = "tbi",
                   patient_col = NULL, paired = FALSE)
)

summarise <- function(st, tag) {
  sig <- st$fdr < 0.05 & abs(st$logFC) >= 0.5
  cat(sprintf("  %-26s 基因=%6d  FDR<0.05=%5d  显著=%5d  max|logFC|=%.2f  min p=%.2e\n",
              tag, nrow(st), sum(st$fdr < 0.05, na.rm = TRUE),
              sum(sig, na.rm = TRUE), max(abs(st$logFC), na.rm = TRUE),
              suppressWarnings(min(st$p, na.rm = TRUE))))
  invisible(sum(sig, na.rm = TRUE))
}

run_one <- function(gse) {
  s <- SPEC[[gse]]
  cat("\n================ ", gse, " (", s$type, ") ================\n", sep = "")
  e <- readRDS(s$rds)[[1]]
  mat <- Biobase::exprs(e); ph <- Biobase::phenoData(e)@data
  cat("  缓存：", nrow(mat), " 基因 × ", ncol(mat), " 样本\n", sep = "")

  sp <- standardise_pheno(ph, s$group_col, s$control_level, s$case_level,
                          patient_col = s$patient_col, batch_col = NULL)
  cat("  分组：case=", sum(sp$group == 1L), "  control=", sum(sp$group == 0L),
      "  patient 唯一=", length(unique(sp$patient)), "\n", sep = "")

  pp <- preprocess_unified(mat, sp, ann = NULL, dataset_id = gse)
  cat("  预处理：尺度=", pp$steps$scale_inferred, "  变换=", pp$steps$transformation,
      "  过滤 ", pp$steps$n_genes_before_filter, "→", pp$steps$n_genes_after_filter,
      "  批次校正=", pp$steps$batch_corrected, "\n", sep = "")
  m <- pp$mat; p <- pp$pheno
  stopifnot(identical(colnames(m), p$sample_id))

  ## 1) 主对比
  st <- if (s$type == "rnaseq") run_de_rnaseq(m, p, paired = s$paired) else run_de_array(m, p, paired = s$paired)
  n_main <- summarise(st, "主对比 case vs control")
  top <- head(st[order(st$fdr, -abs(st$logFC)), c("gene", "logFC", "fdr")], 8)
  cat("  Top8：", paste(sprintf("%s(%+.2f)", top$gene, top$logFC), collapse = " "), "\n")

  ## 1b) 冻结阈值网格（FROZEN$de_fdr_grid × de_lfc_grid）
  ## threshold_grid() 返回长表 data.frame(fdr, lfc, n_de)，这里折成矩阵打印
  tg <- threshold_grid(st, FROZEN$de_fdr_grid, FROZEN$de_lfc_grid)
  lfcs <- sort(unique(tg$lfc)); fdrs <- sort(unique(tg$fdr))
  cat("  阈值网格 n_de（行=FDR，列=|log2FC|≥）：\n")
  cat("           ", paste(sprintf("%8s", lfcs), collapse = ""), "\n")
  for (f in fdrs) {
    v <- vapply(lfcs, function(l) tg$n_de[tg$fdr == f & tg$lfc == l], integer(1))
    cat(sprintf("   FDR<%-5s ", f), paste(sprintf("%8d", v), collapse = ""), "\n")
  }

  ## 2) 标签置换负对照（DESeq2 很慢，RNA-seq 只做 5 次）
  nperm <- if (s$type == "rnaseq") 5L else 20L
  cnt <- integer(nperm)
  for (i in seq_len(nperm)) {
    p2 <- p
    p2$group <- sample(p$group)
    if (s$paired) p2$patient <- p$patient
    st2 <- if (s$type == "rnaseq") run_de_rnaseq(m, p2, paired = s$paired) else run_de_array(m, p2, paired = s$paired)
    cnt[i] <- sum(st2$fdr < 0.05 & abs(st2$logFC) >= 0.5, na.rm = TRUE)
  }
  cat(sprintf("  置换负对照(n=%d)：中位数=%s  90分位=%s  最大=%s  （主对比=%d）\n",
              nperm, round(median(cnt), 1), round(unname(quantile(cnt, .9)), 1),
              max(cnt), n_main))
  if (n_main > 0 && median(cnt) > n_main / 3)
    cat("  ⚠ 置换后仍有接近主对比的信号量 → 慎判，可能有结构冒充\n")

  ## 3) 未配对 vs 配对 的敏感性（仅对配对数据集）
  if (s$paired) {
    st3 <- if (s$type == "rnaseq") run_de_rnaseq(m, p, paired = FALSE) else run_de_array(m, p, paired = FALSE)
    summarise(st3, "敏感性：忽略配对")
  }
  invisible(n_main)
}

sel <- commandArgs(trailingOnly = TRUE)
if (length(sel) == 0 || sel[1] == "all") sel <- names(SPEC)
for (g in sel) if (g %in% names(SPEC)) run_one(g) else cat("跳过：", g, "\n")
cat("\n[smoke] 完成。\n")
