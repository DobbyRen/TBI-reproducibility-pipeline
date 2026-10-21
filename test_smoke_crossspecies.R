## =====================================================================
## test_smoke_crossspecies.R —— 跨物种方案 A + L7 的强信号冒烟测试
##
## 为什么要有这个脚本
##   主线自检（Rscript run_all.R synthetic）用的是两个 human 合成数据集，
##   因此有两条路径【从未被真正执行过】：
##     ① 跨物种基因标识统一（09_ortholog.R）——合成模式直接走 no-op；
##     ② L7 判别性能与锁定操作点——合成数据信号太弱，特征 <3 被设计跳过。
##   本脚本用真实 Ensembl 同源表构造一对"小鼠 / 人"合成数据集，
##   共享全部信号基因，专门把这两条路径跑通并做断言。
##
## 用法：Rscript test_smoke_crossspecies.R
## 退出码：0 = 全部通过；1 = 有断言失败
## =====================================================================

for (f in c("00_config.R", "01_utils.R", "02_io.R", "03_preprocess.R",
            "04_de.R", "05_stability.R", "06_seven_levels.R",
            "07_model_transfer.R", "08_determinants.R", "09_ortholog.R")) {
  source(f, encoding = "UTF-8")
}
ensure_dirs("outputs", "data/raw", "data/processed", "data/gmt")

FAIL <- 0L
chk <- function(label, ok, detail = NULL) {
  d <- if (is.null(detail)) "" else paste(as.character(detail), collapse = ", ")
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", label,
              if (nzchar(d)) paste0(" —— ", d) else ""))
  if (!isTRUE(ok)) FAIL <<- FAIL + 1L
  invisible(ok)
}

## ---------------------------------------------------------------------
## 0. 取真实同源表，构造"小鼠 / 人"共享基因空间
## ---------------------------------------------------------------------
o <- try(fetch_orthologs("mouse", "human"), silent = TRUE)
if (inherits(o, "try-error")) {
  chk("抓取 Ensembl 小鼠→人同源表", FALSE, conditionMessage(attr(o, "condition")))
  quit(status = 1)
}
chk("抓取 Ensembl 小鼠→人同源表", nrow(o$map) > 10000L,
    sprintf("%d 对，Ensembl %s", nrow(o$map), o$release))

set.seed(20261021L)
N_PAIR <- 500L                       # 共享（可映射）基因数
idx <- sample(seq_len(nrow(o$map)), N_PAIR)
mouse_shared <- o$map$from_gene[idx]
human_shared <- o$map$to_gene[idx]

## 填充基因：取未被选中的同源对，保证两侧都是真实 symbol 但互不相同
fill <- setdiff(seq_len(nrow(o$map)), idx)
mouse_fill <- o$map$from_gene[sample(fill, 400)]
human_fill <- o$map$to_gene[sample(fill, 400)]

## ---------------------------------------------------------------------
## 1. 合成两个数据集：小鼠为发现集，人为外部集，信号全部成对共享
## ---------------------------------------------------------------------
## 信号基因必须在两侧取【同一批同源对】，否则重叠只是随机巧合
sig_idx <- sample(seq_len(N_PAIR), 150)
truth_mouse <- mouse_shared[sig_idx]
truth_human <- human_shared[sig_idx]     # 与上一行逐一对应的同源基因

## 效应向量只抽一次、两侧共用 —— 这才是"同一生物学在两个物种中复现"的正确模拟。
## 若两侧各抽各的，效应方向本身就可能相反，L7 迁移失败会是模拟的假象而非代码缺陷。
eff_vec <- rnorm(length(sig_idx), mean = 0, sd = 2.0)
eff_vec[eff_vec > 0] <- eff_vec[eff_vec > 0] + 1.0   # 统一为上调，方向一致
eff_vec[eff_vec < 0] <- eff_vec[eff_vec < 0] - 1.0

make_syn <- function(id, genes_shared, genes_fill, sig_idx, eff_vec, seed, n = 40) {
  set.seed(seed)
  group <- rep(c(0L, 1L), each = n / 2)
  genes <- c(genes_shared, genes_fill)
  p <- length(genes)
  mat <- matrix(rnorm(n * p, mean = 8, sd = 1.2), nrow = p, ncol = n,
                dimnames = list(genes, paste0(id, "_S", seq_len(n))))
  ## 只在共享基因里注入信号，且用较大效应量确保 DE 能检出、L7 特征数 ≥3
  for (j in seq_along(sig_idx)) mat[sig_idx[j], group == 1L] <- mat[sig_idx[j], group == 1L] + eff_vec[j]
  pheno <- data.frame(sample_id = colnames(mat), group = group,
                      patient = paste0(id, "_P", rep(seq_len(n / 2), each = 2)),
                      batch = NA_character_, stringsAsFactors = FALSE)
  list(id = id, mat = mat, pheno = pheno, truth = genes_shared[sig_idx])
}

d_mouse <- make_syn("SYN_MOUSE", mouse_shared, mouse_fill, sig_idx, eff_vec, seed = 101L)
d_human <- make_syn("SYN_HUMAN", human_shared, human_fill, sig_idx, eff_vec, seed = 102L)

## ---------------------------------------------------------------------
## 2. 预处理 → 基因空间统一（方案 A 的关键一步）
## ---------------------------------------------------------------------
pp_m <- preprocess_unified(d_mouse$mat, d_mouse$pheno, ann = NULL, dataset_id = "SYN_MOUSE")
pp_h <- preprocess_unified(d_human$mat, d_human$pheno, ann = NULL, dataset_id = "SYN_HUMAN")
pp_m$data_type <- "array"; pp_h$data_type <- "array"

h_m <- harmonize_gene_space(pp_m$mat, species = "mouse", dataset_id = "SYN_MOUSE")
h_h <- harmonize_gene_space(pp_h$mat, species = "human", dataset_id = "SYN_HUMAN")

chk("小鼠集基因空间统一", identical(h_m$status, "ok"),
    sprintf("%d → %d 行（丢弃 %d，折叠 %d）", h_m$n_in, h_m$n_mapped,
            h_m$n_unmapped_dropped, h_m$n_collapsed))
chk("人集识别为同源空间并原样通过", identical(h_h$status, "native_target_space"),
    sprintf("%d 行", h_h$n_mapped))

pp_m$mat <- h_m$mat; pp_h$mat <- h_h$mat
processed <- list(SYN_MOUSE = pp_m, SYN_HUMAN = pp_h)

chk("统一后共享信号基因可对上", length(intersect(truth_human, rownames(pp_m$mat))) == length(truth_human),
    sprintf("可对上 %d / %d 个真实信号基因", length(intersect(truth_human, rownames(pp_m$mat))),
            length(truth_human)))

## ---------------------------------------------------------------------
## 3. 差异表达 + L4 一致性
## ---------------------------------------------------------------------
st_m <- run_de(pp_m, data_type = "array", paired = FALSE)
st_h <- run_de(pp_h, data_type = "array", paired = FALSE)
set_m <- primary_de(st_m); set_h <- primary_de(st_h)
ov <- length(intersect(set_m, set_h))

chk("两侧均检出 DE 基因", length(set_m) > 20L && length(set_h) > 20L,
    sprintf("小鼠 %d 个 / 人 %d 个", length(set_m), length(set_h)))
## 两侧注入的是同一批同源对，因此重叠应接近满信号数（150）而非随机水平
chk("跨物种 DE 集合重叠 > 0（若未做同源映射此处必为 0）", ov > 80L,
    sprintf("重叠 %d / 信号 %d 个，Jaccard = %.3f", ov, length(truth_mouse),
            jaccard(set_m, set_h)))

l4 <- level4_concordance(st_m, st_h, "SYN_MOUSE", "SYN_HUMAN")
chk("L4 方向一致性可计算", !is.null(l4) && nrow(l4) == 1L,
    if (!is.null(l4)) sprintf("n_overlap = %s", l4$n_overlap[1]) else "")

## ---------------------------------------------------------------------
## 3b. L6 通路层面复制（必须用真实基因符号，合成模式的 G1..G2000 跑不了这步）
## ---------------------------------------------------------------------
h_gmt <- file.path("data/gmt", "msigdb_Hallmark_v7.5.1.gmt")
if (file.exists(h_gmt)) {
  gsets <- read_gmt(h_gmt)
  chk("Hallmark 基因集可读", length(gsets) == 50L, sprintf("%d 个集合", length(gsets)))
  l6 <- try(level6_pathway(st_m, st_h, gsets, "SYN_MOUSE", "SYN_HUMAN", db = "Hallmark"),
            silent = TRUE)
  if (inherits(l6, "try-error")) {
    chk("L6 通路复制可计算", FALSE, conditionMessage(attr(l6, "condition")))
  } else {
    chk("L6 通路复制可计算", !is.null(l6) && nrow(l6) == 1L,
        if (!is.null(l6)) sprintf("公共基因空间 %s，A 显著 %s 个，复制 %s 个，复制率 %s",
                                  l6$n_universe[1], l6$n_sig_A[1], l6$n_replicated[1],
                                  round(l6$replication_rate[1], 3)) else "")
  }
} else {
  message("[跳过] 未找到 ", h_gmt, "，L6 未验证")
}

## ---------------------------------------------------------------------
## 4. L7：判别性能与锁定操作点（此前从未被执行过）
## ---------------------------------------------------------------------
df <- dual_fold_de(pp_m, data_type = "array", paired = FALSE, seed = 20261L)
feat <- intersect(df$internal_set, st_h$gene)
chk("L7 特征数 ≥ 3（满足进入建模的门槛）", length(feat) >= 3L,
    sprintf("内部稳定集 %d 个，与外部集交集 %d 个", length(df$internal_set), length(feat)))

if (length(feat) >= 3L) {
  x_tr <- t(pp_m$mat[feat, , drop = FALSE]); y_tr <- pp_m$pheno$group
  pat  <- pp_m$pheno$patient
  x_te <- t(pp_h$mat[feat, , drop = FALSE]); y_te <- pp_h$pheno$group
  rows <- list()
  for (rule in c(FROZEN$op_rule_primary, FROZEN$op_rule_sensitivity)) {
    md <- try(run_l7(x_tr, y_tr, pat, x_te, y_te, method = FROZEN$model_primary,
                     rule = rule, prevalence = NULL, seed = 20262L,
                     label = "SYN_MOUSE_to_SYN_HUMAN"), silent = TRUE)
    if (inherits(md, "try-error")) {
      chk(paste0("L7 建模（rule = ", rule, "）"), FALSE,
          conditionMessage(attr(md, "condition")))
    } else {
      rows[[length(rows) + 1L]] <- md
      chk(paste0("L7 建模（rule = ", rule, "）"), TRUE,
          sprintf("外部 AUC = %.3f，锁定阈值 = %.3f，外部 sens = %.3f / spec = %.3f，锁定点失效 = %s",
                  md$auc, md$locked_threshold, md$sens, md$spec,
                  md$locked_point_failed))
    }
  }
  if (length(rows)) {
    l7 <- do.call(rbind, rows)
    write_out(l7, "smoke_L7_transfer")
    write_out(summarise_l7(l7), "smoke_L7_summary")
    chk("L7 汇总表可生成", nrow(summarise_l7(l7)) >= 1L)
  }
}

## ---------------------------------------------------------------------
## 5. 收尾
## ---------------------------------------------------------------------
write_out(ortholog_audit_row(h_m), "smoke_ortholog_audit_mouse")
write_out(ortholog_audit_row(h_h), "smoke_ortholog_audit_human")

cat("\n==================== 冒烟测试结果 ====================\n")
cat(if (FAIL == 0L) "全部通过\n" else sprintf("失败 %d 项\n", FAIL))
quit(status = if (FAIL == 0L) 0L else 1L)
