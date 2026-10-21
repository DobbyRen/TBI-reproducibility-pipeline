## =====================================================================
## 02_io.R —— 数据获取、清单读取、合成数据（自检用）
##
## 真实数据：GEOquery 下载 → 缓存为 RDS，避免每次重下。
## 合成数据：RUN_MODE = "synthetic" 时生成两套带部分共享信号的矩阵，
##           用于在联网/GEO 不可用之前先把整条流水线跑通验证。
## =====================================================================

## ---- 读取数据集清单 ----
## 必需列：gse, species, tissue, data_type(array|rnaseq), group_col,
##         control_level, case_level, patient_col, batch_col, paired,
##         time_window, pmid, reported_de_n, claimed_genes
read_manifest <- function(path = "data/datasets_template.csv") {
  if (!file.exists(path)) stop("找不到清单文件：", path)
  m <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("gse", "species", "data_type", "group_col", "control_level", "case_level")
  miss <- setdiff(need, names(m))
  if (length(miss)) stop("清单缺少必需列：", paste(miss, collapse = ", "))
  if (!"paired" %in% names(m)) m$paired <- FALSE
  if (!"patient_col" %in% names(m)) m$patient_col <- NA
  if (!"batch_col" %in% names(m)) m$batch_col <- NA
  m$paired <- toupper(as.character(m$paired)) %in% c("TRUE", "T", "1", "YES")
  m
}

## ---- GEO 下载（带本地缓存）----
fetch_geo <- function(gse, cache_dir = "data/raw") {
  ensure_dirs(cache_dir)
  cache <- file.path(cache_dir, paste0(gse, ".rds"))
  if (file.exists(cache)) {
    message("[io] 命中缓存：", cache)
    return(readRDS(cache))
  }
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("需要 GEOquery：BiocManager::install('GEOquery')")
  }
  esets <- GEOquery::getGEO(gse, destdir = cache_dir, GSEMatrix = TRUE)
  saveRDS(esets, cache)
  esets
}

## ---- 把 ExpressionSet 拆成 矩阵 + 表型 ----
eset_to_parts <- function(eset) {
  list(mat    = Biobase::exprs(eset),
       pheno  = Biobase::phenoData(eset)@data,
       annot  = if (!is.null(Biobase::annotation(eset))) Biobase::annotation(eset) else NA_character_)
}

## ---- 标准化表型：统一列名，避免后续到处猜 ----
## 返回 data.frame：sample_id / group(0=对照,1=损伤) / patient / batch
standardise_pheno <- function(pheno, group_col, control_level, case_level,
                              patient_col = NULL, batch_col = NULL) {
  if (!group_col %in% names(pheno)) stop("表型中没有列：", group_col)
  g_raw <- as.character(pheno[[group_col]])
  ## 用模糊匹配定位水平，避免大小写/空格差异
  hit_ctrl <- grepl(control_level, g_raw, ignore.case = TRUE)
  hit_case <- grepl(case_level,   g_raw, ignore.case = TRUE)
  if (sum(hit_ctrl) == 0 || sum(hit_case) == 0) {
    stop("分组水平匹配失败。现有取值：", paste(unique(g_raw), collapse = " | "))
  }
  group <- ifelse(hit_case, 1L, ifelse(hit_ctrl, 0L, NA_integer_))
  out <- data.frame(
    sample_id = rownames(pheno),
    group     = group,
    stringsAsFactors = FALSE
  )
  out$patient <- if (!is.null(patient_col) && patient_col %in% names(pheno)) {
    as.character(pheno[[patient_col]])
  } else {
    out$sample_id          # 无患者信息时退化为样本水平（需在论文中标注）
  }
  out$batch <- if (!is.null(batch_col) && batch_col %in% names(pheno)) {
    as.character(pheno[[batch_col]])
  } else {
    NA_character_
  }
  out[!is.na(out$group), , drop = FALSE]
}

## ---- 探针注释表：优先平台包，失败则用 ann 数据框 ----
get_annotation <- function(platform, ann_df = NULL) {
  if (!is.null(ann_df)) {
    ann_df <- ann_df[, c("probe_id", "gene_symbol")]
    return(ann_df)
  }
  if (is.na(platform)) return(NULL)
  pkg <- paste0(platform, ".db")
  if (requireNamespace(pkg, quietly = TRUE)) {
    sym <- try(AnnotationDbi::mapIds(get(pkg, asNamespace(pkg)),
                                     keys = AnnotationDbi::keys(get(pkg, asNamespace(pkg))),
                                     column = "SYMBOL", keytype = "PROBEID",
                                     multiVals = "CharacterList"), silent = TRUE)
    if (!inherits(sym, "try-error") && !is.null(sym)) {
      return(data.frame(probe_id = names(sym),
                        gene_symbol = vapply(sym, function(z) {
                          if (length(z) == 0) NA_character_ else paste(z, collapse = "///")
                        }, character(1)),
                        stringsAsFactors = FALSE))
    }
  }
  message("[io] 未找到平台包 ", pkg, "，请自行提供 ann_df（probe_id / gene_symbol）")
  NULL
}

## ---- 缓存 ----
cache_save <- function(obj, name, dir = "data/processed") {
  ensure_dirs(dir)
  f <- file.path(dir, paste0(name, ".rds"))
  saveRDS(obj, f); invisible(f)
}
cache_load <- function(name, dir = "data/processed") {
  f <- file.path(dir, paste0(name, ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

## =====================================================================
## 合成数据：仅供流水线自检，绝不用于论文结果
## =====================================================================
simulate_dataset <- function(id = "SYN_A", n = 40, p = 2000,
                             n_signal = 120, shared_with = NULL,
                             rho_shared = 0.6, seed = 1L) {
  set.seed(seed)
  group <- rep(c(0L, 1L), each = n / 2)
  patient <- paste0(id, "_P", rep(seq_len(n / 2), each = 2))  # 每个患者两份样本
  mat <- matrix(rnorm(n * p, mean = 8, sd = 1.2), nrow = p, ncol = n,
                dimnames = list(paste0("G", seq_len(p)), paste0(id, "_S", seq_len(n))))
  ## 信号基因：与 shared_with 部分共享，用于让 L2/L4 产生非退化数值
  sig <- if (is.null(shared_with)) {
    sample(seq_len(p), n_signal)
  } else {
    k <- round(n_signal * rho_shared)
    c(sample(shared_with, k), sample(setdiff(seq_len(p), shared_with), n_signal - k))
  }
  eff <- rnorm(length(sig), mean = 0, sd = 1.1)
  for (j in seq_along(sig)) {
    mat[sig[j], group == 1L] <- mat[sig[j], group == 1L] + eff[j]
  }
  pheno <- data.frame(sample_id = colnames(mat), group = group,
                      patient = patient, batch = NA_character_,
                      stringsAsFactors = FALSE)
  list(id = id, mat = mat, pheno = pheno, truth_signal = sig)
}
