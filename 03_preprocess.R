## =====================================================================
## 03_preprocess.R —— 统一预处理（对应手稿 §2.5）
##
## 五步，顺序固定，且每一步都写进 audit 表（这张表要作为论文附表逐数据集公布）：
##   (i)   尺度审查      audit_scale()
##   (ii)  注释与探针聚合 collapse_probes()
##   (iii) 变换          log2_if_needed()
##   (iv)  过滤          filter_by_expression()
##   (v)   批次          maybe_combat()
## =====================================================================

## ---- (i) 尺度审查：防止重复取对数或漏取对数 ----
audit_scale <- function(m) {
  v <- as.numeric(m)
  v <- v[is.finite(v)]
  mx <- max(v); mn <- min(v)
  neg <- any(v < 0)
  int_like <- all(abs(v - round(v)) < 1e-8)
  inferred <- if (neg) {
    "log_scale"                       # 有负值 → 已是 log 或比值尺度
  } else if (int_like && mx > 100) {
    "count"                           # 大整数 → 原始计数
  } else if (mx > 50) {
    "linear"                          # 无负值、量级大 → 线性未取对数
  } else {
    "log_scale"
  }
  list(
    inferred     = inferred,
    needs_log2   = inferred %in% c("count", "linear"),
    min = mn, max = mx,
    has_negative = neg,
    integer_like = int_like,
    floor        = unlist(FROZEN$filter_floor[[
      if (inferred == "count") "count" else if (inferred == "linear") "linear" else "log_scale"
    ]])
  )
}

## ---- (ii) 注释与探针聚合 ----
## method = "maxmean": 同基因多探针保留平均表达最高的探针（主分析）
## method = "mean"   : 同基因多探针取均值（敏感性分析）
collapse_probes <- function(mat, ann, method = FROZEN$probe_agg_primary,
                            drop_multi = FROZEN$probe_drop_multi) {
  if (is.null(ann)) return(list(mat = mat, n_dropped_multi = 0L, n_collapsed = 0L))
  ann <- ann[!is.na(ann$gene_symbol) & ann$gene_symbol != "", , drop = FALSE]
  ## 映射到多个基因的探针：直接丢弃
  multi <- grepl("///|;|\\|", ann$gene_symbol)
  n_multi <- sum(multi)
  if (drop_multi && n_multi > 0) ann <- ann[!multi, , drop = FALSE]
  keep <- intersect(rownames(mat), ann$probe_id)
  mat <- mat[keep, , drop = FALSE]
  ann <- ann[match(keep, ann$probe_id), , drop = FALSE]
  sym <- ann$gene_symbol
  dup <- unique(sym[duplicated(sym)])
  if (length(dup) == 0L) {
    rownames(mat) <- sym
    return(list(mat = mat, n_dropped_multi = n_multi, n_collapsed = 0L))
  }
  if (identical(method, "maxmean")) {
    rm_ <- rowMeans(mat, na.rm = TRUE)
    pick <- vapply(dup, function(s) {
      ii <- which(sym == s)
      ii[which.max(rm_[ii])]
    }, integer(1))
    keep_idx <- c(which(!sym %in% dup), as.integer(pick))
  } else {
    mats <- lapply(dup, function(s) {
      ii <- which(sym == s)
      z <- colMeans(mat[ii, , drop = FALSE], na.rm = TRUE)
      z
    })
    m_dup <- do.call(rbind, mats)
    rownames(m_dup) <- dup
    mat_nodup <- mat[!sym %in% dup, , drop = FALSE]
    rownames(mat_nodup) <- sym[!sym %in% dup]
    mat <- rbind(mat_nodup, m_dup)
    return(list(mat = mat, n_dropped_multi = n_multi, n_collapsed = length(dup)))
  }
  mat <- mat[sort(keep_idx), , drop = FALSE]
  rownames(mat) <- sym[sort(keep_idx)]
  list(mat = mat, n_dropped_multi = n_multi, n_collapsed = length(dup))
}

## ---- (iii) 变换：线性/计数 → log2(x + 1) ----
log2_if_needed <- function(mat, needs_log2) {
  if (isTRUE(needs_log2)) {
    if (min(mat, na.rm = TRUE) < 0) stop("线性数据出现负值，无法直接 log2(x+1)，请检查尺度审查结论")
    return(log2(mat + 1))
  }
  mat
}

## ---- (iv) 过滤：至少在一个组中 ≥20% 的样本高于固定下限 ----
filter_by_expression <- function(mat, group, min_frac = FROZEN$filter_min_frac, floor) {
  grp <- unique(group[!is.na(group)])
  keep <- rep(FALSE, nrow(mat))
  for (g in grp) {
    sub <- mat[, group == g, drop = FALSE]
    frac <- rowMeans(sub > floor, na.rm = TRUE)
    keep <- keep | (frac >= min_frac)
  }
  mat[keep, , drop = FALSE]
}

## ---- (v) 批次：仅在显式提供批次变量时才做 ComBat，并保护分组变量 ----
maybe_combat <- function(mat, batch, group, min_per_level = 2L) {
  if (is.null(batch) || all(is.na(batch))) {
    return(list(mat = mat, applied = FALSE, reason = "no_batch_variable_documented"))
  }
  tab <- table(batch)
  if (length(tab) < 2L || min(tab) < min_per_level) {
    return(list(mat = mat, applied = FALSE,
                reason = "batch_has_single_level_or_too_few_samples"))
  }
  if (!requireNamespace("sva", quietly = TRUE)) {
    warning("需要 sva 才能做 ComBat：BiocManager::install('sva')；本次跳过校正")
    return(list(mat = mat, applied = FALSE, reason = "sva_not_installed"))
  }
  mod <- stats::model.matrix(~ group)
  out <- try(sva::ComBat(dat = as.matrix(mat), batch = batch, mod = mod), silent = TRUE)
  if (inherits(out, "try-error")) {
    return(list(mat = mat, applied = FALSE, reason = paste("combat_failed:", conditionMessage(attr(out, "condition")))))
  }
  list(mat = out, applied = TRUE, reason = "documented_batch_variable")
}

## ---- 编排：统一预处理主入口 ----
preprocess_unified <- function(mat, pheno, ann = NULL,
                               probe_method = FROZEN$probe_agg_primary,
                               use_batch = TRUE,
                               dataset_id = "dataset") {
  audit <- audit_scale(mat)
  steps <- list()

  ## (i) 尺度
  steps$scale_inferred <- audit$inferred
  steps$scale_min      <- round(audit$min, 4)
  steps$scale_max      <- round(audit$max, 4)
  steps$needs_log2     <- audit$needs_log2

  ## (ii) 注释与探针聚合
  cp <- collapse_probes(mat, ann, method = probe_method)
  mat <- cp$mat
  steps$n_probes_dropped_multi <- cp$n_dropped_multi
  steps$n_genes_collapsed      <- cp$n_collapsed

  ## (iii) 变换
  mat <- log2_if_needed(mat, audit$needs_log2)
  steps$transformation <- if (audit$needs_log2) "log2(x+1)" else "none"

  ## (iv) 过滤
  n_before <- nrow(mat)
  mat <- filter_by_expression(mat, pheno$group, FROZEN$filter_min_frac, audit$floor)
  steps$filter_floor <- audit$floor
  steps$filter_min_frac <- FROZEN$filter_min_frac
  steps$n_genes_before_filter <- n_before
  steps$n_genes_after_filter  <- nrow(mat)

  ## (v) 批次
  pheno <- pheno[match(colnames(mat), pheno$sample_id), , drop = FALSE]
  if (use_batch) {
    cb <- maybe_combat(mat, pheno$batch, pheno$group)
    mat <- cb$mat
    steps$batch_corrected <- cb$applied
    steps$batch_reason    <- cb$reason
  } else {
    steps$batch_corrected <- FALSE
    steps$batch_reason    <- "disabled_by_argument"
  }

  list(
    dataset_id = dataset_id,
    mat   = mat,
    pheno = pheno,
    audit = audit,
    steps = as.data.frame(steps, stringsAsFactors = FALSE)
  )
}
