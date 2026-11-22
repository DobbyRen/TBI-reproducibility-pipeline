#!/usr/bin/env Rscript
## =====================================================================
## build_geo_cache_generic.R —— GSE89866 / GSE2871 / GSE45997 的缓存构建
##
## 为什么需要它（三个数据集各自的原因，都会写进 §3.7）：
##
##  · GSE89866：RNA-seq，series matrix 无表达值（Sample_type=SRA），
##              数据在 58 个 per-sample HTSeq count 文件里 → 需合并。
##  · GSE2871 / GSE45997：Affymetrix 芯片，series matrix 正常，但矩阵行名是
##              探针 ID。流水线 get_annotation() 拼的包名是 paste0(GPL,".db")
##              （"GPL85.db"/"GPL1355.db"），这类包不存在 → ann=NULL →
##              collapse_probes 空操作 → 后续 L5/L6/跨物种同源映射全线失效。
##              故在【数据准备阶段】用 GEO 平台表把探针折叠成 symbol，
##              聚合规则与 FROZEN$probe_agg_primary / probe_drop_multi 完全一致。
##
## 与 build_geo_cache.R 相同的契约（勿改）：
##   · exprs(eset) 基因 × 样本，列名 = GSM；phenoData 行名相同且顺序一致
##   · annotation = NA_character_（行名已是 symbol，不需再聚合）
##   · 喂【原始尺度】：03 的 audit_scale() 自己判 linear 并做 log2(x+1)
##
## 用法：
##   Rscript build_geo_cache_generic.R --gse GSE89866 [--out data/raw/GSE89866.rds]
##   Rscript build_geo_cache_generic.R --gse GSE2871  --unit HC_4h
##   Rscript build_geo_cache_generic.R --gse GSE45997 --unit ipsi_vs_naive
## =====================================================================
suppressMessages({
  library(Biobase)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}
GSE  <- get_arg("--gse", "GSE89866")
UNIT <- get_arg("--unit", "")
OUT  <- get_arg("--out", file.path("data", "raw", paste0(GSE, ".rds")))
RAWDIR <- "data/raw"

## 预处理脚本目录（脚本在 TBI_repro_pipeline/ 下）
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || !nzchar(here)) here <- "."
if (!file.exists(file.path(here, "00_config.R"))) here <- "."

message("[loader] GSE=", GSE, "  unit=", if (nzchar(UNIT)) UNIT else "(默认)")

## ---------------------------------------------------------------- 公共工具
parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con))
  fields <- list()
  while (TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (!startsWith(line, "!Sample_")) next
    parts <- scan(text = line, what = "", sep = "\t", quote = "\"", quiet = TRUE)
    key <- parts[1]
    fields[[key]][[length(fields[[key]]) + 1L]] <- parts[-1]
  }
  fields
}

## 复刻 03_preprocess.R::collapse_probes（规则一致，勿自创）
collapse_as_pipeline <- function(mat, ann, method = "maxmean", drop_multi = TRUE) {
  ann <- ann[!is.na(ann$gene_symbol) & ann$gene_symbol != "", , drop = FALSE]
  multi <- grepl("///|;|\\|", ann$gene_symbol)
  n_multi <- sum(multi)
  if (drop_multi && n_multi > 0) ann <- ann[!multi, , drop = FALSE]
  keep <- intersect(rownames(mat), ann$probe_id)
  n_nosym <- nrow(mat) - length(keep)
  mat <- mat[keep, , drop = FALSE]
  ann <- ann[match(keep, ann$probe_id), , drop = FALSE]
  sym <- ann$gene_symbol
  dup <- unique(sym[duplicated(sym)])
  n_collapsed <- length(dup)
  if (length(dup) == 0L) {
    rownames(mat) <- sym
    return(list(mat = mat, n_multi = n_multi, n_collapsed = 0L, n_nosym = n_nosym))
  }
  if (identical(method, "maxmean")) {
    rm_ <- rowMeans(mat, na.rm = TRUE)
    pick <- vapply(dup, function(s) {
      ii <- which(sym == s); ii[which.max(rm_[ii])]
    }, integer(1))
    keep_idx <- c(which(!sym %in% dup), as.integer(pick))
    mat <- mat[sort(keep_idx), , drop = FALSE]
    rownames(mat) <- sym[sort(keep_idx)]
  } else {
    stop("本 loader 只实现 maxmean 主分析规则")
  }
  list(mat = mat, n_multi = n_multi, n_collapsed = n_collapsed, n_nosym = n_nosym)
}

finish <- function(mat, pheno, qc_extra = list()) {
  ## 列名 = GSM，pheno 行名 = 同序 GSM（03 依赖这一对齐）
  ord <- order(colnames(mat))
  mat <- mat[, ord, drop = FALSE]
  pheno <- pheno[match(colnames(mat), pheno$gsm), , drop = FALSE]
  rownames(pheno) <- pheno$gsm
  stopifnot(identical(colnames(mat), rownames(pheno)))
  qc <- c(list(gse = GSE, unit = UNIT,
               n_samples = ncol(mat), n_genes = nrow(mat),
               n_case = sum(pheno$group == "tbi"),
               n_control = sum(pheno$group == "control"),
               mat_min = round(min(mat, na.rm = TRUE), 3),
               mat_median = round(median(mat, na.rm = TRUE), 3),
               mat_max = round(max(mat, na.rm = TRUE), 3),
               n_na = sum(is.na(mat)),
               all_integer = all(mat == round(mat)),
               scale_inferred = if (max(mat, na.rm = TRUE) > 50) "linear→log2(x+1) by 03"
                                else "log_scale→none by 03",
               aligned = identical(colnames(mat), rownames(pheno))),
          qc_extra)
  for (k in names(qc)) message(sprintf("  [QC] %-28s %s", k, paste(qc[[k]], collapse = ",")))
  eset <- Biobase::ExpressionSet(assayData = mat,
                                 phenoData = Biobase::AnnotatedDataFrame(pheno),
                                 annotation = NA_character_)
  dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(eset), OUT)
  dir.create("outputs", showWarnings = FALSE)
  tag <- if (nzchar(UNIT)) paste0(GSE, "_", UNIT) else GSE
  write.csv(pheno, file.path("outputs", paste0("qc_samples_", tag, ".csv")), row.names = FALSE)
  write.csv(as.data.frame(qc, stringsAsFactors = FALSE),
            file.path("outputs", paste0("qc_summary_", tag, ".csv")), row.names = FALSE)
  message("[loader] 缓存 → ", OUT, "（", file.size(OUT), " 字节）")
  invisible(NULL)
}

## ================================================================ GSE89866
build_89866 <- function() {
  fmat <- file.path(RAWDIR, "GSE89866_series_matrix.txt.gz")
  fm <- parse_series_matrix(fmat)
  gsm  <- fm$`!Sample_geo_accession`[[1]]
  title <- fm$`!Sample_title`[[1]]
  tp   <- trimws(sub("^time point:", "", fm$`!Sample_characteristics_ch1`[[1]]))
  ## 受试者 = title 去掉末尾 _NN（BL3A01_01 → BL3A01）
  subj <- sub("_[^_]+$", "", title)
  stopifnot(length(gsm) == 58, length(tp) == 58)
  tp_tab <- table(tp)
  if (length(tp_tab) != 2L) stop("时间点不是 2 个：", paste(names(tp_tab), collapse = "|"))
  message("[loader] 时间点：", paste(names(tp_tab), tp_tab, sep = "=", collapse = "  "),
          "；受试者 ", length(unique(subj)), " 位，每人 ",
          paste(unique(table(subj)), collapse = "/"), " 份样本")
  ## 暴露后 = 病例。取值必须互不为子串：base line / 10 days → control / tbi
  group <- ifelse(tp == "10 days", "tbi", ifelse(tp == "base line", "control", NA))
  if (any(is.na(group))) stop("存在未识别的时间点：", paste(unique(tp[is.na(group)]), collapse = "|"))

  d <- file.path(RAWDIR, "GSE89866_counts")
  fs <- list.files(d, pattern = "\\.count\\.txt\\.gz$", full.names = TRUE)
  if (length(fs) != 58) stop("counts 文件数应为 58，实为 ", length(fs))
  key <- sub("_.*$", "", basename(fs))            # GSM2391756
  if (!setequal(key, gsm)) stop("counts 文件名与 GSM 对不上")
  lst <- vector("list", length(fs))
  ref <- NULL
  n_dropped <- 0L
  for (i in seq_along(fs)) {
    dt <- fread(cmd = paste("gzip -dc", shQuote(fs[i])), sep = "\t", header = FALSE,
                stringsAsFactors = FALSE, data.table = FALSE, showProgress = FALSE)
    g <- dt[[1]]; v <- suppressWarnings(as.numeric(dt[[2]]))
    ## HTSeq 末尾的统计行：__too_low_aQual / __not_aligned / __alignment_not_unique
    bad <- startsWith(g, "__")
    n_dropped <- sum(bad)
    g <- g[!bad]; v <- v[!bad]
    if (is.null(ref)) ref <- g else if (!identical(g, ref)) stop("基因顺序不一致：", basename(fs[i]))
    lst[[i]] <- v
  }
  mat <- do.call(cbind, lst)
  rownames(mat) <- ref
  colnames(mat) <- key
  storage.mode(mat) <- "double"
  message("[loader] counts：", nrow(mat), " 基因 × ", ncol(mat),
          " 样本；已剔除 ", n_dropped, " 行 HTSeq 统计行（__ 开头）")

  pheno <- data.frame(gsm = colnames(mat), sample_id = colnames(mat),
                      title = title[match(colnames(mat), gsm)],
                      subject_id = subj[match(colnames(mat), gsm)],
                      timepoint = tp[match(colnames(mat), gsm)],
                      group = group[match(colnames(mat), gsm)],
                      stringsAsFactors = FALSE)
  ## 配对完整性：每位受试者必须恰好一份 control + 一份 tbi
  chk <- tapply(pheno$group, pheno$subject_id, function(z) paste(sort(z), collapse = "+"))
  if (!all(chk == "control+tbi")) stop("配对不完整：", paste(names(chk)[chk != "control+tbi"], collapse = ", "))
  finish(mat, pheno, list(paired_design = TRUE, n_subjects = length(unique(subj)),
                          manifest_patient_col = "subject_id"))
}

## ================================================================ GSE2871
build_2871 <- function() {
  unit <- if (nzchar(UNIT)) UNIT else "HC_4h"
  fmat <- file.path(RAWDIR, "GSE2871_series_matrix.txt.gz")
  fm <- parse_series_matrix(fmat)
  gsm   <- fm$`!Sample_geo_accession`[[1]]
  src   <- fm$`!Sample_source_name_ch1`[[1]]
  desc  <- fm$`!Sample_description`[[1]]
  ## "brain, ipsi Cortex" → side=ipsi / region=Cortex
  side  <- ifelse(grepl("ipsi", src), "ipsi", "contralateral")
  region <- ifelse(grepl("Cortex", src), "Cortex", ifelse(grepl("HC", src), "HC", NA))
  tpt   <- ifelse(grepl("24 hours", desc), "24h", ifelse(grepl("4 hours", desc), "4h", NA))
  samp  <- sub("^Sample(\\d+).*$", "\\1", desc)
  if (any(is.na(region)) || any(is.na(tpt))) stop("无法解析脑区/时间点")

  ## 读表达矩阵
  con <- gzfile(fmat, "rt"); on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  b <- grep("^!series_matrix_table_begin", lines); e <- grep("^!series_matrix_table_end", lines)
  mat <- read.delim(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"),
                    row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(mat); storage.mode(mat) <- "double"
  message("[loader] 探针矩阵：", nrow(mat), " × ", ncol(mat))
  stopifnot(identical(colnames(mat), gsm))

  ## 探针 → symbol
  fannot <- file.path(RAWDIR, "GPL85_symbol.csv")
  if (!file.exists(fannot)) stop("缺 ", fannot, "（先跑 tools/build_platform_annot.py GPL85）")
  ann <- read.csv(fannot, stringsAsFactors = FALSE)
  cp <- collapse_as_pipeline(mat, ann)
  mat <- cp$mat
  message("[loader] 探针折叠：丢弃多映射 ", cp$n_multi, " 个，合并同基因 ",
          cp$n_collapsed, " 组，无 symbol 剔除 ", cp$n_nosym, " 个 → ", nrow(mat), " 基因")

  pheno <- data.frame(gsm = gsm, sample_id = gsm, region = region, side = side,
                      timepoint = tpt, sample_no = samp,
                      group = ifelse(side == "ipsi", "tbi", "control"),
                      stringsAsFactors = FALSE)

  ## 分析单元
  if (unit == "ALL") {
    keep <- rep(TRUE, nrow(pheno))
  } else {
    u <- strsplit(unit, "_")[[1]]
    if (length(u) != 2L) stop("unit 形如 HC_4h / Cortex_24h / ALL，收到：", unit)
    keep <- pheno$region == u[1] & pheno$timepoint == u[2]
  }
  if (sum(keep) == 0) stop("分析单元 ", unit, " 无样本")
  mat <- mat[, keep, drop = FALSE]; pheno <- pheno[keep, , drop = FALSE]
  if (length(unique(pheno$group)) != 2L) stop("分析单元内缺少对照组")

  finish(mat, pheno, list(platform = "GPL85 (RG_U34A)",
                          n_probes_in = cp$n_multi + cp$n_collapsed + nrow(mat) + cp$n_nosym,
                          manifest_patient_col = "NA(样本层，SOFT 未给动物编号)"))
}

## ================================================================ GSE45997
build_45997 <- function() {
  unit <- if (nzchar(UNIT)) UNIT else "ipsi_vs_naive"
  fmat <- file.path(RAWDIR, "GSE45997_series_matrix.txt.gz")
  fm <- parse_series_matrix(fmat)
  gsm  <- fm$`!Sample_geo_accession`[[1]]
  src  <- fm$`!Sample_source_name_ch1`[[1]]
  title <- fm$`!Sample_title`[[1]]
  ch <- fm$`!Sample_characteristics_ch1`
  trt <- NULL
  for (r in ch) if (startsWith(r[1], "treatment:")) trt <- trimws(sub("^treatment:", "", r))
  if (is.null(trt)) stop("找不到 treatment 字段")
  ## 动物编号：题目写法不统一（"Todd 1I" / "2I Todd" / "Xu 13 control"），
  ## 统一抽成 <姓氏><首个数字>：Todd1 / Todd2 / Todd3 / Xu13 / Xu2 / Xu6
  animal <- vapply(title, function(t) {
    fam <- if (grepl("Todd", t)) "Todd" else if (grepl("Xu", t)) "Xu" else ""
    num <- regmatches(t, regexpr("[0-9]+", t))
    paste0(fam, num)
  }, character(1))
  if (any(animal == "")) stop("动物编号解析失败：", paste(title[animal == ""], collapse = ", "))
  message("[loader] 动物编号：", paste(sort(unique(animal)), collapse = ", "))

  con <- gzfile(fmat, "rt"); on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  b <- grep("^!series_matrix_table_begin", lines); e <- grep("^!series_matrix_table_end", lines)
  mat <- read.delim(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"),
                    row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(mat); storage.mode(mat) <- "double"
  message("[loader] 探针矩阵：", nrow(mat), " × ", ncol(mat))
  stopifnot(identical(colnames(mat), gsm))

  fannot <- file.path(RAWDIR, "GPL1355_symbol.csv")
  if (!file.exists(fannot)) stop("缺 ", fannot, "（先跑 tools/build_platform_annot.py GPL1355）")
  ann <- read.csv(fannot, stringsAsFactors = FALSE)
  cp <- collapse_as_pipeline(mat, ann)
  mat <- cp$mat
  message("[loader] 探针折叠：多映射 ", cp$n_multi, "，合并 ", cp$n_collapsed,
          " 组，无 symbol ", cp$n_nosym, " → ", nrow(mat), " 基因")

  cls <- ifelse(trt == "naïve" | trt == "naive", "naive",
         ifelse(grepl("^contralateral", trt), "contralateral", "ipsi"))
  pheno <- data.frame(gsm = gsm, sample_id = gsm, source = src, title = title,
                      treatment = trt, cls = cls, animal = animal,
                      stringsAsFactors = FALSE)
  if (unit == "ipsi_vs_naive") {
    pheno <- pheno[pheno$cls %in% c("ipsi", "naive"), ]
    pheno$group <- ifelse(pheno$cls == "ipsi", "tbi", "control")
    paired <- FALSE
  } else if (unit == "ipsi_vs_contra") {
    pheno <- pheno[pheno$cls %in% c("ipsi", "contralateral"), ]
    pheno$group <- ifelse(pheno$cls == "ipsi", "tbi", "control")
    paired <- TRUE
  } else stop("unit 仅支持 ipsi_vs_naive / ipsi_vs_contra，收到：", unit)
  mat <- mat[, pheno$gsm, drop = FALSE]
  finish(mat, pheno, list(platform = "GPL1355", paired_design = paired,
                          manifest_patient_col = if (paired) "animal" else "NA"))
}

## ---------------------------------------------------------------- 派发
switch(GSE,
  "GSE89866" = build_89866(),
  "GSE2871"  = build_2871(),
  "GSE45997" = build_45997(),
  stop("本脚本未支持：", GSE))
message("[loader] 完成。")
