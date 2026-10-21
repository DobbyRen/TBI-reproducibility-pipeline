## =====================================================================
## run_all.R —— 端到端编排
##
## 两种模式：
##   RUN_MODE = "synthetic"  : 用合成数据把整条流水线跑通（自检用，不产生论文结果）
##   RUN_MODE = "geo"        : 读取 data/datasets_template.csv，真实分析
##
## 用法（在 R 中）：
##   setwd("本文件所在目录")
##   source("run_all.R", encoding = "UTF-8")
##
## 或命令行：
##   Rscript run_all.R synthetic
##   Rscript run_all.R geo
## =====================================================================

## ---- 0. 参数 ----
args <- commandArgs(trailingOnly = TRUE)
RUN_MODE <- if (length(args) >= 1) args[1] else "synthetic"
PREVALENCE_TARGET <- NULL     # 论文需给出目标患病率（如急诊 CT 阳性率），此处留 NULL

## 注意：必须先 source 再建目录——ensure_dirs() 定义在 01_utils.R 里
for (f in c("00_config.R", "01_utils.R", "02_io.R", "03_preprocess.R",
            "04_de.R", "05_stability.R", "06_seven_levels.R",
            "07_model_transfer.R", "08_determinants.R",
            "09_ortholog.R")) {          # 跨物种基因标识统一（方案 A）
  source(f, encoding = "UTF-8")
}

ensure_dirs("outputs", "data/raw", "data/processed", "data/gmt")

log_frozen("00_config.R", "outputs")
message("[run] RUN_MODE = ", RUN_MODE)

## =====================================================================
## 1. 取数：合成 或 GEO
## =====================================================================
if (identical(RUN_MODE, "synthetic")) {
  message("[run] 合成数据模式：仅用于验证流水线可跑通")
  d1 <- simulate_dataset("SYN_A", n = 40, p = 2000, n_signal = 120, seed = 11L)
  d2 <- simulate_dataset("SYN_B", n = 40, p = 2000, n_signal = 120,
                         shared_with = d1$truth_signal, rho_shared = 0.6, seed = 12L)
  manifest <- data.frame(
    gse = c("SYN_A", "SYN_B"), species = c("human", "human"),
    data_type = c("array", "array"), tissue = c("cortex", "blood"),
    n_per_group = c(20, 20), pub_year = c(2021, 2023),
    reported_de_n = c(NA_real_, NA_real_),
    claimed_genes = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  datasets <- list(
    list(id = "SYN_A", mat = d1$mat, pheno = d1$pheno, ann = NULL, data_type = "array"),
    list(id = "SYN_B", mat = d2$mat, pheno = d2$pheno, ann = NULL, data_type = "array")
  )
} else {
  manifest <- read_manifest("data/datasets_template.csv")
  datasets <- lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, ]
    esets <- fetch_geo(row$gse)
    parts <- eset_to_parts(esets[[1]])
    ph <- standardise_pheno(parts$pheno, row$group_col, row$control_level,
                            row$case_level, row$patient_col, row$batch_col)
    list(id = row$gse, mat = parts$mat, pheno = ph,
         ann = get_annotation(parts$annot),
         data_type = if (tolower(row$data_type) %in% c("rnaseq", "rna-seq", "rna_seq"))
                       "rnaseq" else "array")
  })
}

## =====================================================================
## 2. 统一预处理 + 统一差异表达（每个数据集）
## =====================================================================
processed <- list()
for (ds in datasets) {
  message("[run] 预处理：", ds$id)
  pp <- preprocess_unified(ds$mat, ds$pheno, ann = ds$ann, dataset_id = ds$id)
  pp$data_type <- ds$data_type
  processed[[ds$id]] <- pp

  write_out(pp$steps, paste0("audit_", ds$id))
}

## =====================================================================
## 2b. 跨物种基因标识统一（方案 A：Ensembl 一对一同源 → 人基因空间）
##
## 必须放在【差异表达之前】：L4/L5/L6 用的是 DE 基因列表，L5/L7 还要用
## 表达矩阵，若只改列表不改矩阵，两侧会落在不同基因空间里，模型迁移无从谈起。
## 任何统一失败的数据集一律从全部跨数据集层级中剔除——
## 绝不允许"静默退回原基因空间"，那等于把被污染的数字当结果用。
## =====================================================================
species_of <- function(nm) {
  if (!is.null(manifest$species)) {
    s <- manifest$species[match(nm, manifest$gse)]
    if (length(s) == 1 && !is.na(s) && nzchar(as.character(s))) return(as.character(s))
  }
  ORTHO_TARGET          # 清单未写物种时保守地按人处理
}

ortho_rows <- list()
for (nm in names(processed)) {
  pp <- processed[[nm]]
  h  <- harmonize_gene_space(pp$mat, species = species_of(nm),
                             target = ORTHO_TARGET, dataset_id = nm)
  ortho_rows[[nm]] <- ortholog_audit_row(h)
  if (is.null(h$mat)) {
    message("[run][严重] ", nm, " 基因空间统一失败（", h$status,
            "），已从全部跨数据集层级中剔除")
    processed[[nm]] <- NULL
  } else {
    pp$mat <- h$mat
    pp$gene_space <- h
    processed[[nm]] <- pp
    cache_save(pp, paste0("pp_", nm))     # 缓存的是已统一基因空间的版本
  }
}
ortho_tab <- do.call(rbind, ortho_rows)
write_out(ortho_tab, "ortholog_audit")

## 把实际用到的同源映射表导出为 csv，随结果一并交付（透明化 + 可复现）
for (fp in list.files("data/processed", pattern = "^ortholog_.*\\.rds$",
                      full.names = TRUE)) {
  o <- readRDS(fp)
  write_out(data.frame(from_gene = o$map$from_gene, to_gene = o$map$to_gene,
                       ensembl_release = o$release, fetched_at = o$fetched_at,
                       stringsAsFactors = FALSE),
            tools::file_path_sans_ext(basename(fp)))
}
if (any(grepl("^failed_", ortho_tab$status))) {
  message("[run][提醒] 存在基因空间未统一的数据集，详见 outputs/ortholog_audit.csv")
}

stats_list <- list()
for (nm in names(processed)) {
  pp <- processed[[nm]]
  st <- run_de(pp, data_type = pp$data_type, paired = FALSE)
  stats_list[[nm]] <- st
  cache_save(st, paste0("stats_", nm))
}

## ---- L1：DE 负担 ----
l1 <- do.call(rbind, lapply(names(processed), function(nm) {
  rep_n <- suppressWarnings(as.numeric(manifest$reported_de_n[manifest$gse == nm]))
  level1_burden(stats_list[[nm]], nrow(processed[[nm]]$mat),
                reported_de_n = if (length(rep_n) == 0) NA_real_ else rep_n[1],
                dataset_id = nm)
}))
write_out(l1, "L1_de_burden")

l1_grid <- do.call(rbind, lapply(names(stats_list), function(nm) {
  level1_threshold_burden(stats_list[[nm]], nm)
}))
write_out(l1_grid, "L1_threshold_grid")

## ---- 内部双层拆分（§2.7）与负对照 ----
stability <- list()
for (nm in names(processed)) {
  pp <- processed[[nm]]
  message("[run] 双层拆分：", nm)
  df <- dual_fold_de(pp, data_type = pp$data_type, paired = FALSE,
                     seed = seed_for(paste0("split_", nm)))
  stability[[nm]] <- df
  cache_save(df, paste0("stability_", nm))
}

stab_tab <- do.call(rbind, lapply(names(stability), function(nm) {
  s <- stability[[nm]]
  data.frame(dataset_id = nm,
             n_full_de = s$n_full,
             n_fold1 = s$n_fold1, n_fold2 = s$n_fold2,
             n_intersect = s$n_intersect, n_concordant = s$n_concordant,
             prop_survived = s$prop_survived,
             stringsAsFactors = FALSE)
}))
write_out(stab_tab, "internal_stability")

## 负对照：规模匹配随机抽样（每个数据集互为外部）
if (length(stats_list) >= 2L) {
  nms <- names(stats_list)
  neg <- list()
  for (a in nms) {
    for (b in setdiff(nms, a)) {
      universe <- stats_list[[a]]$gene
      rd <- size_matched_draw(universe, length(stability[[a]]$internal_set),
                              seed = seed_for(paste0("neg_", a, "_", b)))
      neg[[paste0(a, "_vs_", b)]] <- negctrl_compare(
        stability[[a]]$internal_set, rd, stats_list[[b]],
        seed = seed_for(paste0("negcmp_", a, "_", b)))
    }
  }
  neg_tab <- do.call(rbind, lapply(names(neg), function(k) {
    d <- neg[[k]]; d$pair <- k; d
  }))
  write_out(neg_tab, "negative_control")
}

## =====================================================================
## 3. L2–L4：成对比较（需 ≥2 个数据集）
## =====================================================================
if (length(stats_list) >= 2L) {
  nms <- names(stats_list)
  universe <- Reduce(union, lapply(stats_list, function(s) s$gene))
  sets <- lapply(stats_list, primary_de)

  l2 <- level2_matrix(sets, universe)
  write_out(l2, "L2_membership")

  l3a <- do.call(rbind, lapply(seq_along(nms), function(i) {
    do.call(rbind, lapply(seq_along(nms), function(j) {
      if (j <= i) return(NULL)
      level3_robustness(stats_list[[nms[i]]], stats_list[[nms[j]]], nms[i], nms[j])
    }))
  }))
  write_out(l3a, "L3_threshold_robustness")

  l3b <- do.call(rbind, lapply(seq_along(nms), function(i) {
    do.call(rbind, lapply(seq_along(nms), function(j) {
      if (j <= i) return(NULL)
      level3_topk(stats_list[[nms[i]]], stats_list[[nms[j]]], labelA = nms[i], labelB = nms[j])
    }))
  }))
  write_out(l3b, "L3_topk_overlap")

  l4 <- do.call(rbind, lapply(seq_along(nms), function(i) {
    do.call(rbind, lapply(seq_along(nms), function(j) {
      if (j <= i) return(NULL)
      level4_concordance(stats_list[[nms[i]]], stats_list[[nms[j]]], nms[i], nms[j])
    }))
  }))
  write_out(l4, "L4_concordance")
}

## =====================================================================
## 4. L5：预先指定基因的确认率
## =====================================================================
l5 <- list()
for (a in names(stats_list)) {
  for (b in setdiff(names(stats_list), a)) {
    claimed <- NULL
    if ("claimed_genes" %in% names(manifest)) {
      raw <- manifest$claimed_genes[manifest$gse == a]
      if (length(raw) == 1 && !is.na(raw) && nzchar(raw)) {
        claimed <- trimws(strsplit(as.character(raw), ";")[[1]])
      }
    }
    if (is.null(claimed)) {
      ## 合成/占位：用发现集 top 30 作为"声明基因"，仅用于跑通流程
      claimed <- head(rank_genes(stats_list[[a]]), 30L)
    }
    l5[[paste0(a, "_claimed_in_", b)]] <-
      level5_stratified(claimed, stats_list[[b]], study_id = a)
  }
}
l5_tab <- do.call(rbind, lapply(names(l5), function(k) {
  d <- l5[[k]]; d$pair <- k; d
}))
write_out(l5_tab, "L5_confirmation")

## =====================================================================
## 5. L6：通路层面复制（需要 .gmt 文件；放在 data/gmt/ 下）
## =====================================================================
gmt_files <- list.files("data/gmt", pattern = "\\.gmt$", full.names = TRUE)
if (length(gmt_files) > 0L && length(stats_list) >= 2L) {
  nms <- names(stats_list)
  l6 <- list()
  for (f in gmt_files) {
    gs <- read_gmt(f)
    db <- tools::file_path_sans_ext(basename(f))
    for (i in seq_along(nms)) for (j in seq_along(nms)) {
      if (j <= i) next
      r <- level6_pathway(stats_list[[nms[i]]], stats_list[[nms[j]]], gs,
                          nms[i], nms[j], db = db)
      if (!is.null(r)) l6[[length(l6) + 1L]] <- r
    }
  }
  if (length(l6)) write_out(do.call(rbind, l6), "L6_pathway")
} else {
  message("[run] 未发现 data/gmt/*.gmt，跳过 L6。",
          "请从 MSigDB 下载 Hallmark / KEGG / GO_BP 的 .gmt 放入该目录。")
}

## =====================================================================
## 6. L7：判别性能与锁定操作点迁移
## =====================================================================
if (length(processed) >= 2L) {
  nms <- names(processed)
  disc <- nms[1]; ext <- nms[2]
  feat <- intersect(stability[[disc]]$internal_set, stats_list[[ext]]$gene)
  message("[run] L7：发现集 ", disc, " → 外部集 ", ext, "；特征数 = ", length(feat))

  if (length(feat) >= 3L) {
    x_tr <- t(processed[[disc]]$mat[feat, , drop = FALSE])
    y_tr <- processed[[disc]]$pheno$group
    pat  <- processed[[disc]]$pheno$patient
    x_te <- t(processed[[ext]]$mat[feat, , drop = FALSE])
    y_te <- processed[[ext]]$pheno$group

    rows <- list()
    for (rule in c(FROZEN$op_rule_primary, FROZEN$op_rule_sensitivity)) {
      md <- try(run_l7(x_tr, y_tr, pat, x_te, y_te,
                       method = FROZEN$model_primary, rule = rule,
                       prevalence = PREVALENCE_TARGET,
                       seed = seed_for(paste0("L7_", rule)),
                       label = paste0(disc, "_to_", ext)), silent = TRUE)
      if (inherits(md, "try-error")) {
        message("[run] L7 主模型失败：", conditionMessage(attr(md, "condition")))
      } else {
        rows[[length(rows) + 1L]] <- md
      }
    }
    if (length(rows)) {
      l7 <- do.call(rbind, rows)
      write_out(l7, "L7_transfer")
      write_out(summarise_l7(l7), "L7_summary")
      message("[run] 锁定操作点失效比例 = ",
              round(mean(l7$locked_point_failed), 3))
    }
  } else {
    message("[run] 两数据集公共特征不足 3 个，跳过 L7")
  }
}

## =====================================================================
## 7. 影响因素（真实数据模式下才有意义）
## =====================================================================
if (!identical(RUN_MODE, "synthetic")) {
  X <- build_predictors(manifest)
  ok <- stats::complete.cases(X)
  if (sum(ok) >= 8L) {
    ## 结局必须与预测变量【同长度、同顺序】。
    ## L1 是一研究一行，可直接对齐；L5 目前是一"研究×外部数据集×分层"多行，
    ## 需先按研究聚合（例如取各外部数据集确认率的中位数）再进入模型。
    l5_by_study <- stats::aggregate(
      confirm_rate ~ study_id,
      data = l5_tab[l5_tab$stratum == "all_claimed", , drop = FALSE],
      FUN = function(z) stats::median(z, na.rm = TRUE))
    outcomes <- list(
      L1_prop_de = l1$prop_de[ok],
      L5_confirm_rate = l5_by_study$confirm_rate[
        match(manifest$gse[ok], l5_by_study$study_id)]
    )
    det <- run_determinants(outcomes, X[ok, , drop = FALSE],
                            cluster = manifest$gse[ok])
    if (!is.null(det)) write_out(det, "determinants")
  } else {
    message("[run] 完整预测变量不足 8 个研究，跳过影响因素分析")
  }
}

## =====================================================================
## 8. 收尾
## =====================================================================
writeLines(capture.output(utils::sessionInfo()), "outputs/session_info.txt")
message("[run] 完成。输出目录：outputs/")
message("[run] 提醒：L5 的方向判据目前固定为 sign(logFC)==+1，",
        "请按各数据集的实际对比方向统一校准后再出正式结果。")
