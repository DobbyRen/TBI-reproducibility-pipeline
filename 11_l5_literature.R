## =====================================================================
## 11_l5_literature.R —— H5 真实版：文献结论级点名基因在独立数据集中的确认率
##
## 与 10-27 代理版的区别（这是本脚本存在的唯一理由）：
##   10-27 的 L5_confirmation.csv 用的是「各数据集自身 top-30 DEG」当作被声明基因
##   （run_all.R:258 的占位分支），那是同源自证，不是对文献声明的外部验证。
##   本脚本改用双人提取 + 仲裁定版的文献结论级点名基因（data/claimed_genes_literature.csv）。
##
## 独立性：每个研究只在其【未使用过】的数据集中验证（exclude_datasets；
##         S021 与 GSE104687 同属 Allen/ACT 队列，同源即非独立，已排除）。
## 判定：目标数据集中 FDR < 0.05 且 |logFC| >= 0.5 才算「确认」（主判据，冻结参数）。
## =====================================================================

for (f in c("00_config.R", "01_utils.R", "02_io.R", "03_preprocess.R",
            "04_de.R", "05_stability.R", "06_seven_levels.R",
            "07_model_transfer.R", "08_determinants.R", "09_ortholog.R")) {
  source(f, encoding = "UTF-8")
}
ensure_dirs("outputs", "data/raw", "data/processed", "data/gmt")

claim <- read.csv("data/claimed_genes_literature.csv",
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
stats <- read.csv("outputs/stats_all_genes.csv",
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")

## ---- 同源映射表（原物种 symbol → 人 symbol），Ensembl 一对一同源 ----
map_m <- read.csv("outputs/ortholog_mmusculus_to_hsapiens.csv",
                  stringsAsFactors = FALSE)
map_r <- read.csv("outputs/ortholog_rnorvegicus_to_hsapiens.csv",
                  stringsAsFactors = FALSE)
mm <- setNames(map_m$to_gene, map_m$from_gene)
rr <- setNames(map_r$to_gene, map_r$from_gene)

## 文献基因 → 人基因空间
##   1) 人：直接大写
##   2) 鼠/大鼠：先查同源表；未命中则回退「首字母大写→全大写」并标记为弱映射
to_human <- function(g, sp) {
  sp <- tolower(sp)
  if (sp %in% c("human", "homo sapiens")) return(list(gene = toupper(g), how = "direct"))
  tbl <- if (grepl("mouse|mus", sp)) mm else if (grepl("rat|rattus", sp)) rr else NULL
  if (!is.null(tbl) && g %in% names(tbl)) return(list(gene = tbl[[g]], how = "ortholog"))
  if (g %in% names(tbl)) return(list(gene = tbl[[g]], how = "ortholog"))
  ## 同源表未命中：大小写回退（弱映射，单独标记，敏感性分析中可剔除）
  return(list(gene = toupper(g), how = "case_fallback"))
}

claim$h_gene <- NA_character_
claim$map_how <- NA_character_
for (i in seq_len(nrow(claim))) {
  r <- to_human(trimws(claim$gene_symbol_raw[i]), claim$species[i])
  claim$h_gene[i] <- r$gene
  claim$map_how[i] <- r$how
}

## ---- 目标数据集（有 DE 结果的 6 个）----
targets <- sort(unique(stats$dataset_id))

## 每个数据集的显著基因集合（主判据）
sig_of <- lapply(targets, function(d) {
  s <- stats[stats$dataset_id == d, ]
  s$sigid <- !is.na(s$fdr) & s$fdr < FROZEN$de_fdr_primary &
    abs(s$logFC) >= FROZEN$de_lfc_primary
  s
})
names(sig_of) <- targets

## ---- 逐「研究 × 目标数据集」验证 ----
res <- list()
for (sid in unique(claim$study_id)) {
  sub <- claim[claim$study_id == sid, ]
  excl <- trimws(strsplit(sub$exclude_datasets[1], ";")[[1]])
  excl <- excl[nzchar(excl)]
  for (tier in c("primary", "secondary")) {
    g_all <- sub$h_gene[sub$tier == tier & !is.na(sub$h_gene)]
    g_all <- unique(g_all)
    if (length(g_all) == 0) next
    for (d in setdiff(targets, excl)) {
      s <- sig_of[[d]]
      idx <- match(g_all, s$gene)
      hit <- !is.na(idx)
      n_map <- sum(hit)
      confirmed <- sum(hit & s$sigid[idx[hit]])
      lfc <- s$logFC[idx[hit]]
      res[[length(res) + 1L]] <- data.frame(
        study_id = sid, tier = tier, target_dataset = d,
        n_claimed = length(g_all), n_in_target_space = n_map,
        n_confirmed = confirmed,
        confirm_rate = if (n_map > 0) confirmed / n_map else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
}
res <- do.call(rbind, res)
write.csv(res, "outputs/L5_literature_confirmation.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

## ---- 汇总（primary 为主报告，secondary 仅敏感性）----
ci_wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  p <- x / n
  z <- qnorm(0.975)
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - hw), min(1, ctr + hw))
}

p <- res[res$tier == "primary", ]
tot_x <- sum(p$n_confirmed, na.rm = TRUE)
tot_n <- sum(p$n_in_target_space, na.rm = TRUE)
ci <- ci_wilson(tot_x, tot_n)
message(sprintf("\n[H5-汇总] 文献点名基因 %d 条研究；可映射「研究×独立数据集」单元 %d 个",
                length(unique(p$study_id)), nrow(p)))
message(sprintf("[H5-汇总] 总体确认率 = %d / %d = %.4f  (95%%CI %.4f–%.4f)",
                tot_x, tot_n, tot_x / tot_n, ci[1], ci[2]))

## 按目标数据集分层（关键：0-DEG 数据集会把总体拉低，必须单独呈现）
by_ds <- do.call(rbind, lapply(targets, function(d) {
  q <- p[p$target_dataset == d, ]
  if (nrow(q) == 0) return(NULL)
  x <- sum(q$n_confirmed, na.rm = TRUE); n <- sum(q$n_in_target_space, na.rm = TRUE)
  cc <- ci_wilson(x, n)
  data.frame(target_dataset = d, n_studies = nrow(q), n_genes_mapped = n,
             n_confirmed = x, confirm_rate = if (n > 0) x / n else NA_real_,
             ci_lo = cc[1], ci_hi = cc[2], stringsAsFactors = FALSE)
}))
write.csv(by_ds, "outputs/L5_literature_by_dataset.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

## 按研究汇总（每条研究在其全部独立数据集上的合并确认率）
by_st <- do.call(rbind, lapply(sort(unique(p$study_id)), function(sid) {
  q <- p[p$study_id == sid, ]
  x <- sum(q$n_confirmed, na.rm = TRUE); n <- sum(q$n_in_target_space, na.rm = TRUE)
  cc <- ci_wilson(x, n)
  data.frame(study_id = sid, n_targets = nrow(q),
             n_genes_mapped = n, n_confirmed = x,
             confirm_rate = if (n > 0) x / n else NA_real_,
             ci_lo = cc[1], ci_hi = cc[2], stringsAsFactors = FALSE)
}))
write.csv(by_st, "outputs/L5_literature_by_study.csv", row.names = FALSE,
          fileEncoding = "UTF-8")

## 排除 0-DEG 数据集后的敏感性总率（GSE104687 与 GSE160763 主对比下 0 个 DEG）
nz <- targets[sapply(targets, function(d) sum(sig_of[[d]]$sigid, na.rm = TRUE) > 0)]
q <- p[p$target_dataset %in% nz, ]
x2 <- sum(q$n_confirmed, na.rm = TRUE); n2 <- sum(q$n_in_target_space, na.rm = TRUE)
cc2 <- ci_wilson(x2, n2)
message(sprintf("[H5-敏感性] 仅在有 DEG 产出的数据集(%s)上：%d/%d = %.4f (95%%CI %.4f–%.4f)",
                paste(nz, collapse = ","), x2, n2, x2 / n2, cc2[1], cc2[2]))

## 映射方式统计
message("\n[H5] 基因映射方式分布：")
print(table(claim$map_how[claim$tier == "primary"]))

message("\n[H5] 按目标数据集：")
print(by_ds)
message("\n[H5] 按研究：")
print(by_st)
