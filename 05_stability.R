## =====================================================================
## 05_stability.R —— 内部特征稳定性（对应手稿 §2.7，REDE-2Fold 的 TBI 适配）
##
## 三条纪律：
##  1. 拆分单元是【患者/动物】，不是样本。同一患者的样本必须整体进入同一折。
##  2. 两折各自独立做 DE，只保留"两折都达到主阈值 且 方向一致"的基因。
##  3. 必须做负对照：规模匹配的随机抽样。否则无法回答"双层拆分到底有没有用"
##     —— 审稿人一定会问。
## =====================================================================

## ---- 患者水平分层拆分 ----
patient_level_split <- function(pheno, prop = FROZEN$split_prop,
                                stratify = FROZEN$split_stratify, seed) {
  set.seed(seed)
  pats <- unique(pheno$patient)
  ## 每个患者的分组（取该患者样本的众数组）
  pat_grp <- vapply(pats, function(pp) {
    v <- pheno$group[pheno$patient == pp]
    as.integer(names(sort(table(v), decreasing = TRUE))[1])
  }, integer(1))

  f1 <- character(0)
  if (stratify) {
    for (g in unique(pat_grp)) {
      pp <- pats[pat_grp == g]
      k <- max(1L, floor(length(pp) * prop))
      f1 <- c(f1, sample(pp, k))
    }
  } else {
    f1 <- sample(pats, floor(length(pats) * prop))
  }
  f2 <- setdiff(pats, f1)

  i1 <- which(pheno$patient %in% f1)
  i2 <- which(pheno$patient %in% f2)

  ## 断言：绝不出现同一患者跨折
  stopifnot(length(intersect(pheno$patient[i1], pheno$patient[i2])) == 0L)
  list(
    fold1 = i1, fold2 = i2,
    n_patient_f1 = length(f1), n_patient_f2 = length(f2),
    n_sample_f1 = length(i1), n_sample_f2 = length(i2)
  )
}

## ---- 双层拆分差异表达 ----
dual_fold_de <- function(obj, data_type = "array", paired = FALSE, seed) {
  sp <- patient_level_split(obj$pheno, seed = seed)
  res1 <- run_de(list(mat = obj$mat[, sp$fold1, drop = FALSE],
                      pheno = obj$pheno[sp$fold1, , drop = FALSE]),
                 data_type = data_type, paired = paired)
  res2 <- run_de(list(mat = obj$mat[, sp$fold2, drop = FALSE],
                      pheno = obj$pheno[sp$fold2, , drop = FALSE]),
                 data_type = data_type, paired = paired)

  full <- run_de(obj, data_type = data_type, paired = paired)
  full_set <- primary_de(full)

  s1 <- primary_de(res1); s2 <- primary_de(res2)
  both <- intersect(s1, s2)

  ## 方向一致：这是"内部稳定"的定义，缺了它就只是交集
  d1 <- sign(res1$logFC[match(both, res1$gene)])
  d2 <- sign(res2$logFC[match(both, res2$gene)])
  keep <- both[d1 == d2 & d1 != 0]

  list(
    internal_set = keep,
    n_fold1 = length(s1), n_fold2 = length(s2),
    n_intersect = length(both), n_concordant = length(keep),
    n_full = length(full_set),
    prop_survived = if (length(full_set) > 0) length(keep) / length(full_set) else NA_real_,
    split = sp,
    fold1_stats = res1, fold2_stats = res2, full_stats = full,
    full_set = full_set
  )
}

## ---- 负对照：规模匹配的随机抽样 ----
## 从全队列排名中随机抽 n 个基因（不做任何拆分/筛选），用于回答
## "双层拆分的收益是来自拆分策略，还是仅仅来自集合变小？"
size_matched_draw <- function(universe, n, seed) {
  set.seed(seed)
  if (length(universe) == 0L || n <= 0L) return(character(0))
  sample(universe, min(n, length(universe)))
}

## ---- 负对照比较：两组的外部确认率之差 ----
## internal_genes / random_genes 各自在 external_stats 中的确认率
## （确认 = 达到主阈值 且 方向一致，与 L5 定义一致）
negctrl_compare <- function(internal_genes, random_genes, external_stats, seed) {
  confirm <- function(genes, ext) {
    ii <- match(genes, ext$gene)
    ii <- ii[!is.na(ii)]
    if (length(ii) == 0L) return(rep(0L, 0))
    sig <- ext$fdr[ii] < FROZEN$de_fdr_primary & abs(ext$logFC[ii]) >= FROZEN$de_lfc_primary
    as.integer(if (FROZEN$require_direction) (sig & sign(ext$logFC[ii]) == 1L) else sig)
  }
  ci <- confirm(internal_genes, external_stats)
  cr <- confirm(random_genes, external_stats)

  ## 用"目标方向为上调"作为方向判据；若原数据集以对照为基线，则下调基因为负值。
  ## 论文中需按各数据集的实际对比方向统一取反，此处演示用 +1。
  b_int <- boot_ci_prop(ci, seed = seed)
  b_rnd <- boot_ci_prop(cr, seed = seed)
  d     <- boot_ci_diff(ci, cr, seed = seed)

  data.frame(
    n_internal = length(internal_genes), n_random = length(random_genes),
    confirm_rate_internal = b_int[["est"]],
    ci_lo_internal = b_int[["lo"]], ci_hi_internal = b_int[["hi"]],
    confirm_rate_random = b_rnd[["est"]],
    ci_lo_random = b_rnd[["lo"]], ci_hi_random = b_rnd[["hi"]],
    diff = d[["diff"]], diff_lo = d[["lo"]], diff_hi = d[["hi"]],
    stringsAsFactors = FALSE
  )
}

## ---- 多次重抽的负对照分布（稳健性：不只跑一次）----
negctrl_repeated <- function(internal_genes, universe, external_stats,
                             draws = FROZEN$negctrl_draws, seed) {
  set.seed(seed)
  n <- length(internal_genes)
  rates <- numeric(draws)
  for (b in seq_len(draws)) {
    rg <- size_matched_draw(universe, n, seed = seed_for(paste0("neg", b)))
    ii <- match(rg, external_stats$gene)
    ii <- ii[!is.na(ii)]
    if (length(ii) == 0L) { rates[b] <- NA_real_; next }
    sig <- external_stats$fdr[ii] < FROZEN$de_fdr_primary &
           abs(external_stats$logFC[ii]) >= FROZEN$de_lfc_primary
    if (FROZEN$require_direction) sig <- sig & sign(external_stats$logFC[ii]) == 1L
    rates[b] <- mean(sig)
  }
  c(mean = mean(rates, na.rm = TRUE),
    lo = as.numeric(quantile(rates, 0.025, na.rm = TRUE)),
    hi = as.numeric(quantile(rates, 0.975, na.rm = TRUE)))
}
