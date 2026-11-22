#!/usr/bin/env Rscript
## =====================================================================
## build_geo_cache_sc.R —— 两个单细胞数据集（GSE269748 / GSE160763）的缓存构建
##
## 为什么单列一个脚本：
##   · 这两个是 10X 单细胞，GEO 只提供 per-cell mtx，没有 bulk 表达矩阵；
##     counts 由 tools/pseudobulk.py 按【样本】聚合生成（预先指定 2026-10-27）。
##   · 流水线一行不改：仍走 fetch_geo() 的 data/raw/<GSE>.rds 缓存通道，
##     下游 03 判 count 尺度 → run_de_rnaseq (DESeq2)。
##
## 预先指定（早于任何单细胞 DE 结果，已写入 data/datasets_template.csv notes）：
##   · GSE269748：主单元 = 24h CCI ipsilateral (9) vs 24h Naive (10)
##                 contralateral 不进主对比（同鼠对侧半球）
##                 7d / 6m 无时间匹配对照 → 不做 DE，只描述
##   · GSE160763：2×2 设计，injury Control(4) vs TBI(4)；n<5 → 降级敏感性
##
## 用法：
##   Rscript build_geo_cache_sc.R --gse GSE269748 --unit 24h_ipsi_vs_naive
##   Rscript build_geo_cache_sc.R --gse GSE269748 --unit 24h_male_ipsi_vs_naive
##   Rscript build_geo_cache_sc.R --gse GSE160763 --unit injury_main
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
GSE  <- get_arg("--gse", "GSE269748")
UNIT <- get_arg("--unit", "")
PB   <- get_arg("--pb", file.path("data", "raw", paste0("pb_", GSE, ".csv")))
OUT  <- get_arg("--out", file.path("data", "raw", paste0(GSE, ".rds")))
RAWDIR <- "data/raw"

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0 || !nzchar(here)) here <- "."
if (!file.exists(file.path(here, "00_config.R"))) here <- "."

message("[sc-loader] GSE=", GSE, "  unit=", if (nzchar(UNIT)) UNIT else "(默认)")

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

finish <- function(mat, pheno, qc_extra = list()) {
  ## 与 build_geo_cache_generic.R 完全相同的契约
  ord <- order(colnames(mat))
  mat <- mat[, ord, drop = FALSE]
  pheno <- pheno[match(colnames(mat), pheno$gsm), , drop = FALSE]
  rownames(pheno) <- pheno$gsm
  stopifnot(identical(colnames(mat), rownames(pheno)))
  stopifnot(!any(is.na(pheno$group)))
  stopifnot(all(pheno$group %in% c("control", "tbi")))   ## 互不为子串，standardise_pheno 用 grepl
  qc <- c(list(gse = GSE, unit = UNIT,
               n_samples = ncol(mat), n_genes = nrow(mat),
               n_case = sum(pheno$group == "tbi"),
               n_control = sum(pheno$group == "control"),
               mat_min = round(min(mat, na.rm = TRUE), 3),
               mat_median = round(median(mat, na.rm = TRUE), 3),
               mat_max = round(max(mat, na.rm = TRUE), 3),
               n_na = sum(is.na(mat)),
               all_integer = all(mat == round(mat)),
               scale_inferred = "count→DESeq2 by 03",
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
  message("[sc-loader] 缓存 → ", OUT, "（", file.size(OUT), " 字节）")
  invisible(NULL)
}

read_pb <- function() {
  if (!file.exists(PB)) stop("缺 pseudobulk 输出 ", PB, "（先跑 tools/pseudobulk.py）")
  dt <- fread(PB, sep = ",", header = TRUE, stringsAsFactors = FALSE,
              data.table = FALSE, showProgress = FALSE)
  g <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  rownames(mat) <- g
  storage.mode(mat) <- "double"
  ## 列名形如 GSM8326331_RNASEQ01 → 只保留 GSM
  colnames(mat) <- sub("_.*$", "", colnames(mat))
  if (anyDuplicated(colnames(mat))) stop("pseudobulk 列名去重后出现重复 GSM")
  message("[sc-loader] pseudobulk：", nrow(mat), " 基因 × ", ncol(mat), " 样本")
  if (!all(mat == round(mat))) stop("pseudobulk 应为整数 counts")
  mat
}

## ================================================================ GSE269748
build_269748 <- function() {
  unit <- if (nzchar(UNIT)) UNIT else "24h_ipsi_vs_naive"
  fmat <- file.path(RAWDIR, "GSE269748_series_matrix.txt.gz")
  fm <- parse_series_matrix(fmat)
  gsm  <- fm$`!Sample_geo_accession`[[1]]
  title <- fm$`!Sample_title`[[1]]
  ch <- fm$`!Sample_characteristics_ch1`
  trt <- NULL
  for (r in ch) if (startsWith(r[1], "treatment:")) trt <- trimws(sub("^treatment:", "", r))
  if (is.null(trt)) stop("找不到 treatment 字段")

  ## 时间点：title 前缀 24 hours / 7 day / 6 month
  tp <- ifelse(grepl("^24 hours", title), "24h",
        ifelse(grepl("^7 day", title), "7d",
        ifelse(grepl("^6 month", title), "6m", NA)))
  if (any(is.na(tp))) stop("时间点解析失败：", paste(title[is.na(tp)], collapse = " | "))
  ## 性别：⚠️ "female" 含 "male"，必须先判 female
  sex <- ifelse(grepl("female", title), "F", ifelse(grepl("male", title), "M", NA))
  ## 侧别：⚠️ 原始 title 有拼写错误 "CCI ipisilateral replicate 4"（ipisilateral），
  ##         grepl("ipsi") 命中不了它 → 改用「明确 contralateral 才是对侧」的判定
  side <- ifelse(grepl("contralateral", title), "contralateral", "ipsilateral")
  arm <- ifelse(grepl("^CCI$", trt), "CCI",
         ifelse(grepl("^Naive$", trt), "Naive",
         ifelse(grepl("CCI *\\+ *HS|Combined", trt), "CCI_HS",
         ifelse(grepl("Repetitive|weight drop", trt), "rCHI", NA))))
  if (any(is.na(arm))) stop("处理臂解析失败：", paste(unique(trt[is.na(arm)]), collapse = " | "))

  mat <- read_pb()
  if (!setequal(colnames(mat), gsm)) {
    miss <- setdiff(gsm, colnames(mat)); extra <- setdiff(colnames(mat), gsm)
    stop("pseudobulk 样本与 series matrix 对不上：缺 ", paste(miss, collapse = ","),
         "；多 ", paste(extra, collapse = ","))
  }
  pheno <- data.frame(gsm = colnames(mat),
                      title = title[match(colnames(mat), gsm)],
                      timepoint = tp[match(colnames(mat), gsm)],
                      arm = arm[match(colnames(mat), gsm)],
                      side = side[match(colnames(mat), gsm)],
                      sex = sex[match(colnames(mat), gsm)],
                      stringsAsFactors = FALSE)

  ## 分析单元（预先指定）
  group <- rep(NA_character_, nrow(pheno))
  if (unit == "24h_ipsi_vs_naive") {
    keep <- pheno$timepoint == "24h" &
            ((pheno$arm == "CCI" & pheno$side == "ipsilateral") | pheno$arm == "Naive")
    group[keep] <- ifelse(pheno$arm[keep] == "CCI", "tbi", "control")
    paired <- FALSE
  } else if (unit == "24h_male_ipsi_vs_naive") {      ## 敏感性①
    keep <- pheno$timepoint == "24h" & pheno$sex == "M" &
            ((pheno$arm == "CCI" & pheno$side == "ipsilateral") | pheno$arm == "Naive")
    group[keep] <- ifelse(pheno$arm[keep] == "CCI", "tbi", "control")
    paired <- FALSE
  } else if (unit == "24h_ipsi_vs_contra") {          ## 敏感性②：同鼠配对，n=3
    keep <- pheno$timepoint == "24h" & pheno$arm == "CCI" & pheno$side %in%
            c("ipsilateral", "contralateral")
    ## 只保留有配对对象的鼠：ipsilateral replicate 4/5/6 ↔ contralateral 1/2/3
    keep <- keep & (pheno$side == "contralateral" |
                    grepl("replicate *(4|5|6)", pheno$title))
    group[keep] <- ifelse(pheno$side[keep] == "ipsilateral", "tbi", "control")
    paired <- TRUE
  } else if (unit == "24h_CCIHS_vs_naive") {          ## 敏感性③
    keep <- pheno$timepoint == "24h" & pheno$arm %in% c("CCI_HS", "Naive")
    group[keep] <- ifelse(pheno$arm[keep] == "CCI_HS", "tbi", "control")
    paired <- FALSE
  } else if (unit == "24h_rCHI_vs_naive") {           ## 敏感性④
    keep <- pheno$timepoint == "24h" & pheno$arm %in% c("rCHI", "Naive")
    group[keep] <- ifelse(pheno$arm[keep] == "rCHI", "tbi", "control")
    paired <- FALSE
  } else stop("GSE269748 不支持的 unit：", unit)

  pheno$group <- group
  mat <- mat[, !is.na(group), drop = FALSE]
  pheno <- pheno[!is.na(group), , drop = FALSE]
  if (sum(pheno$group == "tbi") < 2 || sum(pheno$group == "control") < 2)
    stop("分析单元样本不足：tbi=", sum(pheno$group == "tbi"),
         " control=", sum(pheno$group == "control"))
  if (paired) {
    ## 配对号：series matrix 写明 contralateral replicate n 与 ipsilateral replicate n+3 同鼠
    ## （"same mouse as ipsilateral replicate 5"）→ 统一折算成鼠号 1/2/3
    ## ⚠️ 陷阱：对侧样本的 title 里有两个 "replicate"
    ##    （"contralateral replicate 2  (same mouse as ipsilateral replicate 5)"），
    ##    用贪婪的 sub(".*replicate *([0-9]+).*") 会取到后一个 → 必须用 regexpr 取第一个
    m <- regmatches(pheno$title, regexpr("replicate *[0-9]+", pheno$title))
    rep_no <- as.integer(sub("replicate *", "", m))
    if (any(is.na(rep_no))) stop("配对号解析失败：", paste(pheno$title[is.na(rep_no)], collapse = " | "))
    mouse <- ifelse(pheno$side == "ipsilateral", rep_no - 3L, rep_no)
    if (anyDuplicated(paste(pheno$side, mouse))) stop("配对号重复")
    chk <- tapply(pheno$side, mouse, function(z) paste(sort(z), collapse = "+"))
    if (!all(chk == "contralateral+ipsilateral"))
      stop("配对不完整：", paste(names(chk)[chk != "contralateral+ipsilateral"], collapse = ", "))
    pheno$subject_id <- paste0("mouse", mouse)
    message("[sc-loader] 配对：", paste(names(chk), collapse = ", "), " 各 ipsi+contra")
  }
  finish(mat, pheno, list(paired_design = paired,
                          manifest_patient_col = if (paired) "subject_id" else "NA",
                          note = "10X pseudobulk by sample"))
}

## ================================================================ GSE160763
build_160763 <- function() {
  unit <- if (nzchar(UNIT)) UNIT else "injury_main"
  fmat <- file.path(RAWDIR, "GSE160763_series_matrix.txt.gz")
  fm <- parse_series_matrix(fmat)
  gsm  <- fm$`!Sample_geo_accession`[[1]]
  title <- fm$`!Sample_title`[[1]]
  ch <- fm$`!Sample_characteristics_ch1`
  getc <- function(pref) {
    for (r in ch) if (startsWith(r[1], pref)) return(trimws(sub(paste0("^", pref), "", r)))
    NULL
  }
  inj <- getc("injury:"); trt <- getc("treatment:"); sex <- getc("Sex:")
  tis <- getc("tissue:")
  if (is.null(inj) || is.null(trt)) stop("找不到 injury / treatment 字段")

  mat <- read_pb()
  if (!setequal(colnames(mat), gsm)) stop("pseudobulk 样本与 series matrix 对不上")
  pheno <- data.frame(gsm = colnames(mat),
                      title = title[match(colnames(mat), gsm)],
                      injury = inj[match(colnames(mat), gsm)],
                      treatment = trt[match(colnames(mat), gsm)],
                      sex = if (is.null(sex)) NA_character_ else sex[match(colnames(mat), gsm)],
                      tissue = if (is.null(tis)) NA_character_ else tis[match(colnames(mat), gsm)],
                      stringsAsFactors = FALSE)

  ## 分析单元
  if (unit == "injury_main") {
    ## 2×2 全部 8 样本；主效应 injury，treatment 作协变量（DESeq2 设计 ~treatment + injury）
    pheno$group <- ifelse(pheno$injury == "TBI", "tbi",
                   ifelse(pheno$injury == "Control", "control", NA))
    paired <- FALSE
  } else if (unit == "vehicle_only") {
    ## 敏感性：仅 Vehicle 亚组 2 vs 2（n<5，仅参考不作推断）
    keep <- pheno$treatment == "Vehicle"
    pheno <- pheno[keep, , drop = FALSE]
    pheno$group <- ifelse(pheno$injury == "TBI", "tbi", "control")
    mat <- mat[, pheno$gsm, drop = FALSE]
    paired <- FALSE
  } else stop("GSE160763 不支持的 unit：", unit)

  if (any(is.na(pheno$group))) stop("分组解析失败：", paste(unique(pheno$injury[is.na(pheno$group)]), collapse = "|"))
  finish(mat, pheno, list(paired_design = paired,
                          manifest_patient_col = "NA",
                          covariate = if (unit == "injury_main") "treatment" else "NA",
                          note = "10X pseudobulk by sample; n=4/4 -> 敏感性"))
}

## ---------------------------------------------------------------- 派发
switch(GSE,
  "GSE269748" = build_269748(),
  "GSE160763" = build_160763(),
  stop("本脚本未支持：", GSE))
message("[sc-loader] 完成。")
