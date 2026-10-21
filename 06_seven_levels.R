## =====================================================================
## 06_seven_levels.R —— 七级可重复性度量 L1–L6（对应手稿 §2.8）
## L7（判别性能与锁定操作点迁移）在 07_model_transfer.R
##
## 证据层级规则：每一级只在上一级成立时才解读下一级。
##   只有 L5 成立才可以用 "replicated gene" 这个表述。
##   L6 必须区分"通路复制"与"通路内基因复制"——两者含义不同。
## =====================================================================

## ================= L1 预处理有效性与差异基因负担 =================
## 分子：统一阈值下的 DE 基因数；分母：该数据集通过过滤的基因总数
level1_burden <- function(stats, n_filtered, reported_de_n = NA_real_,
                          dataset_id = "dataset") {
  n_de <- length(primary_de(stats))
  data.frame(
    dataset_id   = dataset_id,
    n_genes_tested = n_filtered,
    n_de         = n_de,
    prop_de      = if (n_filtered > 0) n_de / n_filtered else NA_real_,
    reported_de_n = reported_de_n,
    ratio_recomputed_over_reported = if (!is.na(reported_de_n) && reported_de_n > 0) {
      n_de / reported_de_n
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}

## L1 附带：阈值网格上的 DE 负担（原始研究阈值不同会怎样）
level1_threshold_burden <- function(stats, dataset_id = "dataset") {
  tab <- threshold_grid(stats)
  tab$dataset_id <- dataset_id
  tab[, c("dataset_id", "fdr", "lfc", "n_de")]
}

## ================= L2 清单精确重叠 =================
level2_membership <- function(setA, setB, universe, labelA = "A", labelB = "B") {
  j <- jaccard(setA, setB)
  oc <- overlap_coef(setA, setB)
  p <- hyper_p(setA, setB, universe)
  data.frame(
    pair = paste0(labelA, "_vs_", labelB),
    n_A = length(unique(setA)), n_B = length(unique(setB)),
    n_intersect = length(intersect(setA, setB)),
    jaccard = j, overlap_coef = oc, hyper_p = p,
    stringsAsFactors = FALSE
  )
}

## pairwise 矩阵（Fig. 2 用 UpSet / 热图呈现集合结构）
level2_matrix <- function(sets, universe) {
  nm <- names(sets)
  out <- NULL
  for (i in seq_along(nm)) {
    for (j in seq_along(nm)) {
      if (j <= i) next
      out <- rbind(out, level2_membership(sets[[i]], sets[[j]], universe, nm[i], nm[j]))
    }
  }
  out
}

## ================= L3 阈值稳健性与高排名重叠 =================
## 在固定阈值网格上重算，考察结论是否随"合理但不同"的阈值翻转
level3_robustness <- function(statsA, statsB, labelA = "A", labelB = "B") {
  ga <- threshold_grid(statsA); gb <- threshold_grid(statsB)
  sa <- attr(ga, "sets"); sb <- attr(gb, "sets")
  keys <- intersect(names(sa), names(sb))
  universe <- union(statsA$gene, statsB$gene)
  rows <- lapply(keys, function(k) {
    data.frame(threshold = k,
               n_A = length(sa[[k]]), n_B = length(sb[[k]]),
               jaccard = jaccard(sa[[k]], sb[[k]]),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab$pair <- paste0(labelA, "_vs_", labelB)
  tab[, c("pair", "threshold", "n_A", "n_B", "jaccard")]
}

## top-k 交集：k = 25 / 50 / 100
level3_topk <- function(statsA, statsB, ks = FROZEN$topk_grid,
                        labelA = "A", labelB = "B") {
  ra <- rank_genes(statsA); rb <- rank_genes(statsB)
  data.frame(
    pair = paste0(labelA, "_vs_", labelB),
    k = ks,
    overlap = vapply(ks, function(k) topk_overlap(ra, rb, k), numeric(1)),
    stringsAsFactors = FALSE
  )
}

## ================= L4 效应方向一致性 =================
level4_concordance <- function(statsA, statsB, labelA = "A", labelB = "B") {
  common <- intersect(statsA$gene, statsB$gene)
  ia <- match(common, statsA$gene); ib <- match(common, statsB$gene)

  ## 全基因组 log2FC 的 Spearman ρ
  rho <- suppressWarnings(stats::cor(statsA$logFC[ia], statsB$logFC[ib],
                                     method = "spearman"))

  ## 两数据集中均显著的基因
  sigA <- primary_de(statsA); sigB <- primary_de(statsB)
  both <- intersect(sigA, sigB)
  if (length(both) == 0L) {
    return(data.frame(pair = paste0(labelA, "_vs_", labelB),
                      n_common = length(common), spearman_rho = rho,
                      n_double_sig = 0L, concordant = NA_integer_,
                      concordance_rate = NA_real_, sign_test_p = NA_real_,
                      stringsAsFactors = FALSE))
  }
  da <- sign(statsA$logFC[match(both, statsA$gene)])
  db <- sign(statsB$logFC[match(both, statsB$gene)])
  conc <- sum(da == db)
  bt <- stats::binom.test(conc, length(both), p = 0.5,
                          alternative = FROZEN$sign_test_alt)
  data.frame(
    pair = paste0(labelA, "_vs_", labelB),
    n_common = length(common), spearman_rho = rho,
    n_double_sig = length(both), concordant = conc,
    concordance_rate = conc / length(both),
    sign_test_p = bt$p.value,
    stringsAsFactors = FALSE
  )
}

## ================= L5 预先指定基因的确认率 =================
## 分子：在独立数据集中达到主阈值【且方向一致】的预先指定基因数
## 分母：该研究声明的基因家族总数（多重性校正的分母即为该家族）
## 必须分层报告：全部声明基因 vs |log2FC| >= 1 的大效应子集
level5_confirmation <- function(claimed_genes, external_stats,
                                lfc_min = NULL, study_id = "study") {
  ii <- match(unique(claimed_genes), external_stats$gene)
  found <- !is.na(ii)
  ii <- ii[found]
  g <- external_stats$gene[ii]
  lfc <- external_stats$logFC[ii]
  fdr <- external_stats$fdr[ii]

  sig <- fdr < FROZEN$de_fdr_primary & abs(lfc) >= FROZEN$de_lfc_primary
  if (!is.null(lfc_min)) sig <- sig & abs(lfc) >= lfc_min
  if (FROZEN$require_direction) sig <- sig & sign(lfc) == 1L   # 见下方方向说明

  b <- boot_ci_prop(as.integer(sig), seed = seed_for(paste0("L5", study_id)))
  data.frame(
    study_id = study_id,
    n_claimed = length(unique(claimed_genes)),
    n_mapped  = sum(found),
    n_confirmed = sum(sig),
    confirm_rate = b[["est"]],
    ci_lo = b[["lo"]], ci_hi = b[["hi"]],
    lfc_min = if (is.null(lfc_min)) NA_real_ else lfc_min,
    stringsAsFactors = FALSE
  )
}

## 分层封装：全部 + 大效应子集，两个数都要给
level5_stratified <- function(claimed_genes, external_stats, study_id = "study") {
  rbind(
    cbind(stratum = "all_claimed",
          level5_confirmation(claimed_genes, external_stats, NULL, study_id)),
    cbind(stratum = paste0("abs_logFC_ge_", FROZEN$l5_lfc_large),
          level5_confirmation(claimed_genes, external_stats,
                              FROZEN$l5_lfc_large, study_id))
  )
}

## ================= L6 通路与程序层面的复制 =================
## 基因集：优先 clusterProfiler，退化为自实现的超几何 ORA（数学等价）
read_gmt <- function(path) {
  if (!file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  lapply(strsplit(ln, "\t"), function(z) {
    if (length(z) < 3) return(NULL)
    list(name = z[1], genes = unique(z[3:length(z)]))
  }) -> lst
  lst <- Filter(Negate(is.null), lst)
  out <- lapply(lst, function(z) z$genes)
  names(out) <- vapply(lst, function(z) z$name, character(1))
  out
}

ora_hyper <- function(genes, universe, gsets,
                      min_set = FROZEN$enrich_min_set,
                      max_set = FROZEN$enrich_max_set) {
  universe <- unique(universe)
  genes <- unique(genes)
  gsets <- gsets[vapply(gsets, length, integer(1)) >= min_set]
  gsets <- gsets[vapply(gsets, length, integer(1)) <= max_set]
  if (length(gsets) == 0L) return(NULL)
  n_uni <- length(universe); k <- length(genes)
  res <- lapply(names(gsets), function(nm) {
    s <- intersect(gsets[[nm]], universe)
    m <- length(s)
    q <- length(intersect(genes, s))
    p <- if (m == 0 || k == 0) NA_real_ else
         stats::phyper(q - 1, m, n_uni - m, k, lower.tail = FALSE)
    data.frame(term = nm, set_size = m, overlap = q, p = p,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, res)
  tab$fdr <- stats::p.adjust(tab$p, method = "BH")
  tab[order(tab$p), ]
}

## 富集方向：以重叠基因的平均 logFC 符号定义（ORA 没有 NES，方向需显式给出）
enrich_direction <- function(genes, gsets, stats) {
  lfc <- stats$logFC[match(genes, stats$gene)]
  vapply(names(gsets), function(nm) {
    ii <- intersect(gsets[[nm]], genes)
    if (length(ii) == 0L) return(NA_real_)
    mean(lfc[match(ii, genes)], na.rm = TRUE)
  }, numeric(1))
}

level6_pathway <- function(statsA, statsB, gsets,
                           labelA = "A", labelB = "B",
                           db = "Hallmark") {
  if (is.null(gsets) || length(gsets) == 0L) {
    message("[L6] 没有可用的基因集（请把 .gmt 放到 data/gmt/）")
    return(NULL)
  }
  ## 基因空间必须显式记录：两数据集的公共基因空间
  universe <- intersect(statsA$gene, statsB$gene)
  setA <- primary_de(statsA); setB <- primary_de(statsB)

  eA <- ora_hyper(setA, universe, gsets)
  eB <- ora_hyper(setB, universe, gsets)
  if (is.null(eA) || is.null(eB)) return(NULL)

  ## 方向必须在各自数据集内计算，不能共用同一个 stats
  dA <- enrich_direction(setA, gsets, statsA)
  dB <- enrich_direction(setB, gsets, statsB)

  eA$direction <- dA[match(eA$term, names(dA))]
  eB$direction <- dB[match(eB$term, names(dB))]

  ## fdr 可能含 NA（如某通路在公共基因空间里为空、超几何无法计算），
  ## 直接比较会得到 NA，进而让后面的 sum() 返回 NA、if (NA > 0) 报错。
  ## 这里显式把 NA 判为非显著，并用 na.rm 求和。
  eA$sig <- !is.na(eA$fdr) & eA$fdr < FROZEN$enrich_fdr
  eB$sig <- !is.na(eB$fdr) & eB$fdr < FROZEN$enrich_fdr

  rep_tab <- merge(eA, eB, by = "term", suffixes = c("_A", "_B"))
  same_dir <- !is.na(rep_tab$direction_A) & !is.na(rep_tab$direction_B) &
    (sign(rep_tab$direction_A) == sign(rep_tab$direction_B))
  rep_tab$replicated <- rep_tab$sig_A & rep_tab$sig_B & same_dir

  ## 关键区分：通路复制 ≠ 通路内基因复制
  rep_tab$gene_jaccard <- vapply(rep_tab$term, function(nm) {
    jaccard(intersect(setA, gsets[[nm]]), intersect(setB, gsets[[nm]]))
  }, numeric(1))

  n_sig_a <- sum(eA$sig, na.rm = TRUE)
  n_sig_b <- sum(eB$sig, na.rm = TRUE)
  n_rep   <- sum(rep_tab$replicated, na.rm = TRUE)

  data.frame(
    database = db,
    pair = paste0(labelA, "_vs_", labelB),
    n_universe = length(universe),
    n_sig_A = n_sig_a, n_sig_B = n_sig_b,
    n_replicated = n_rep,
    replication_rate = if (n_sig_a > 0) n_rep / n_sig_a else NA_real_,
    median_gene_jaccard = median(rep_tab$gene_jaccard[rep_tab$replicated], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
