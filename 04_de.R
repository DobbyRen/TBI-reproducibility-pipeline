## =====================================================================
## 04_de.R —— 统一差异表达（对应手稿 §2.6 (i)）
##
## 微阵列 → limma；RNA-seq 计数 → DESeq2（失败/缺失时退化为 limma-voom）。
## 两套方法都返回同一张标准表：gene / logFC / AveExpr / p / fdr
## 后续所有层级只认这张表 —— 这是"统一重算"能成立的前提。
## =====================================================================

## ---- 微阵列：limma ----
run_de_array <- function(mat, pheno, paired = FALSE) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("需要 limma：BiocManager::install('limma')")
  y <- as.integer(pheno$group)
  if (paired && !is.null(pheno$patient)) {
    design <- stats::model.matrix(~ factor(pheno$patient) + y)
  } else {
    design <- stats::model.matrix(~ y)
  }
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = "y", sort.by = "none", number = Inf,
                        adjust.method = "BH")
  out <- data.frame(gene = rownames(tt), logFC = tt$logFC,
                    AveExpr = tt$AveExpr, p = tt$P.Value, fdr = tt$adj.P.Val,
                    stringsAsFactors = FALSE)
  out[order(out$fdr, -abs(out$logFC), out$gene), ]
}

## ---- RNA-seq：DESeq2 ----
run_de_rnaseq <- function(counts, pheno, paired = FALSE) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    message("[de] DESeq2 不可用，退化为 limma-voom")
    return(run_de_voom(counts, pheno, paired = paired))
  }
  counts <- round(as.matrix(counts))
  grp <- factor(ifelse(pheno$group == 1L, "case", "ctrl"), levels = c("ctrl", "case"))
  cdat <- data.frame(row.names = colnames(counts), group = grp)
  if (paired && !is.null(pheno$patient)) {
    cdat$patient <- factor(pheno$patient)
    design <- ~ patient + group
  } else {
    design <- ~ group
  }
  dds <- try(DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                            colData = cdat, design = design), silent = TRUE)
  if (inherits(dds, "try-error")) return(run_de_voom(counts, pheno, paired = paired))
  dds <- try(DESeq2::DESeq(dds, quiet = TRUE), silent = TRUE)
  if (inherits(dds, "try-error")) return(run_de_voom(counts, pheno, paired = paired))
  res <- DESeq2::results(dds, contrast = c("group", "case", "ctrl"))
  res <- as.data.frame(res)
  out <- data.frame(gene = rownames(res), logFC = res$log2FoldChange,
                    AveExpr = res$baseMean, p = res$pvalue, fdr = res$padj,
                    stringsAsFactors = FALSE)
  out$fdr[is.na(out$fdr)] <- 1
  out[order(out$fdr, -abs(out$logFC), out$gene), ]
}

## ---- RNA-seq 兜底：limma-voom ----
run_de_voom <- function(counts, pheno, paired = FALSE) {
  if (!requireNamespace("limma", quietly = TRUE)) stop("需要 limma：BiocManager::install('limma')")
  y <- as.integer(pheno$group)
  design <- if (paired && !is.null(pheno$patient)) {
    stats::model.matrix(~ factor(pheno$patient) + y)
  } else {
    stats::model.matrix(~ y)
  }
  v <- limma::voom(as.matrix(counts), design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = "y", sort.by = "none", number = Inf, adjust.method = "BH")
  out <- data.frame(gene = rownames(tt), logFC = tt$logFC,
                    AveExpr = tt$AveExpr, p = tt$P.Value, fdr = tt$adj.P.Val,
                    stringsAsFactors = FALSE)
  out[order(out$fdr, -abs(out$logFC), out$gene), ]
}

## ---- 统一入口：按数据类型自动分派 ----
run_de <- function(obj, data_type = c("array", "rnaseq"), paired = FALSE) {
  data_type <- match.arg(data_type)
  if (data_type == "array") run_de_array(obj$mat, obj$pheno, paired = paired)
  else                      run_de_rnaseq(obj$mat, obj$pheno, paired = paired)
}

## ---- 阈值网格：L1 的 DE 负担 + L3 的稳健性 ----
## 返回长表：每个 (fdr, lfc) 组合下的 DE 基因数与基因集合（集合放在 attr 里）
threshold_grid <- function(stats, fdr_grid = FROZEN$de_fdr_grid,
                           lfc_grid = FROZEN$de_lfc_grid) {
  rows <- list()
  sets <- list()
  for (f in fdr_grid) {
    for (l in lfc_grid) {
      g <- stats$gene[stats$fdr < f & abs(stats$logFC) >= l]
      rows[[length(rows) + 1L]] <- data.frame(fdr = f, lfc = l, n_de = length(g),
                                              stringsAsFactors = FALSE)
      sets[[paste0("fdr", f, "_lfc", l)]] <- g
    }
  }
  tab <- do.call(rbind, rows)
  attr(tab, "sets") <- sets
  tab
}

## ---- WGCNA 模块（固定软阈值规则）—— L6 与文章 II 会用到 ----
rebuild_modules <- function(mat, rsq_cut = 0.85, power_cap = 20L, min_size = 30L) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    message("[wgcna] 未安装，跳过模块重建")
    return(NULL)
  }
  powers <- 1:power_cap
  sft <- WGCNA::pickSoftThreshold(t(mat), powerVector = powers, verbose = 0)
  ok <- which(sft$fitIndices[, "SFT.R.sq"] >= rsq_cut)
  power <- if (length(ok)) min(powers[ok]) else power_cap
  net <- WGCNA::blockwiseModules(t(mat), power = power, minModuleSize = min_size,
                                 verbose = 0, numericLabels = FALSE)
  list(power = power, modules = net$colors, sft = sft$fitIndices)
}
