#!/usr/bin/env Rscript
## =====================================================================
## build_geo_cache.R —— 为"GEO 无 processed series matrix"的数据集造缓存
##
## 背景：run_all.R 的 geo 分支走 fetch_geo() → GEOquery::getGEO()。
## GSE104687 / GSE89866 的 series matrix 里没有表达值（Sample_type=SRA），
## getGEO 拿不到矩阵。但 fetch_geo() 的第一件事就是查 data/raw/<GSE>.rds
## 缓存——所以这里直接把「矩阵 + 表型」装成 ExpressionSet 存成缓存，
## 流水线代码一个字都不用改。
##
## 用法：
##   Rscript tools/build_geo_cache.R --gse GSE104687 [--variant rin|tbt] \
##          [--region HIP|PCx|FWM|TCx|ALL] [--out data/raw/GSE104687.rds]
##
## 契约（与 02_io.R / 03_preprocess.R 对齐，勿改）：
##   · exprs(eset)  基因 × 样本，列名 = 样本 ID（用 GSM）
##   · phenoData    行名 = 同样的 GSM，且【顺序与矩阵列一致】
##   · annotation   = NA  → get_annotation() 返回 NULL → collapse_probes 不聚合
##                    （本数据集行名已是 gene symbol，本就不需要探针聚合）
##   · 矩阵喂【原始 FPKM】：03 的 audit_scale() 会自己判 linear 并做 log2(x+1)
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

GSE     <- get_arg("--gse", "GSE104687")
VARIANT <- get_arg("--variant", "rin")        # rin = RIN_corrected（主分析）, tbt = TBT_normalized（敏感性）
REGION  <- get_arg("--region", "HIP")         # 主分析单元；ALL = 四区合并
OUT     <- get_arg("--out", file.path("data", "raw", paste0(GSE, ".rds")))
RAWDIR  <- "data/raw"

message("[loader] GSE=", GSE, "  variant=", VARIANT, "  region=", REGION)

## ---------------------------------------------------------------- 1. 元数据
parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
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

pick_ch <- function(fields, prefix) {
  rows <- fields[["!Sample_characteristics_ch1"]]
  if (is.null(rows)) return(NULL)
  for (r in rows) {
    if (length(r) && startsWith(r[1], paste0(prefix, ":"))) {
      return(trimws(sub(paste0("^", prefix, ":"), "", r)))
    }
  }
  NULL
}

fmat  <- file.path(RAWDIR, paste0(GSE, "_series_matrix.txt.gz"))
if (!file.exists(fmat)) stop("找不到 series matrix：", fmat)
fm <- parse_series_matrix(fmat)

gsm     <- fm$`!Sample_geo_accession`[[1]]
title   <- fm$`!Sample_title`[[1]]
desc    <- fm$`!Sample_description`[[1]]      # = FPKM 列头的数字型样本 ID
region  <- fm$`!Sample_source_name_ch1`[[1]]
donor   <- pick_ch(fm, "donor_name")
batch   <- pick_ch(fm, "batch")
tbi     <- pick_ch(fm, "ever_tbi_w_loc")
dem     <- pick_ch(fm, "act_demented")
sex     <- pick_ch(fm, "Sex")
age     <- pick_ch(fm, "age_at_death")
rin     <- pick_ch(fm, "rin")
braak   <- pick_ch(fm, "braak")
cerad   <- pick_ch(fm, "cerad")
apo_e4  <- pick_ch(fm, "apo_e4_allele")
ntbi    <- pick_ch(fm, "num_tbi_w_loc")

n <- length(gsm)
stopifnot(length(desc) == n, length(region) == n, length(tbi) == n)
message("[loader] 元数据样本数 = ", n, "；供者 = ", length(unique(donor)),
        "；脑区 = ", paste(sort(unique(region)), collapse = "/"))

## 分组列取值必须互不为子串：02_io.R 用 grepl 模糊匹配，
## "no_tbi" 与 "tbi" 会同时命中 → 全被判成病例。故用 control / tbi。
tbi_group <- ifelse(toupper(trimws(tbi)) %in% c("Y", "YES", "TRUE"), "tbi", "control")
message("[loader] 分组：tbi = ", sum(tbi_group == "tbi"), "，control = ",
        sum(tbi_group == "control"))

## ---------------------------------------------------------------- 2. 表达矩阵
fexpr <- switch(VARIANT,
  rin = file.path(RAWDIR, paste0(GSE, "_RIN_corrected_FPKM.csv.gz")),
  tbt = file.path(RAWDIR, paste0(GSE, "_TBT_normalized_FPKM.csv.gz")),
  stop("未知 variant：", VARIANT))
if (!file.exists(fexpr)) stop("找不到表达文件：", fexpr)
message("[loader] 读取 ", basename(fexpr), "（56 MB 级，需十几秒）")

dt <- tryCatch(
  fread(cmd = paste("gzip -dc", shQuote(fexpr)), sep = ",", header = TRUE,
        stringsAsFactors = FALSE, data.table = FALSE, showProgress = FALSE),
  error = function(e) {
    message("[loader] fread 失败，回退 read.csv：", conditionMessage(e))
    read.csv(gzfile(fexpr), stringsAsFactors = FALSE, check.names = FALSE)
  })

genes <- dt[[1]]
cols  <- colnames(dt)[-1]
mat   <- as.matrix(dt[, -1, drop = FALSE])
storage.mode(mat) <- "double"
rownames(mat) <- genes
message("[loader] 矩阵：", nrow(mat), " 基因 × ", ncol(mat), " 样本")

## 列头（数字 ID）→ GSM
if (!setequal(cols, desc)) {
  stop("FPKM 列头与 !Sample_description 对不上：",
       "仅列头有 ", length(setdiff(cols, desc)), " 个，仅元数据有 ",
       length(setdiff(desc, cols)), " 个")
}
colnames(mat) <- gsm[match(cols, desc)]
if (anyDuplicated(rownames(mat))) {
  dupn <- table(rownames(mat))
  message("[loader] 行名有 ", sum(dupn > 1), " 个重复 symbol，按均值合并")
  agg  <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  mat  <- agg / as.numeric(dupn[rownames(agg)])     # 按行名匹配的计数相除 = 均值
}
mat <- mat[, order(colnames(mat)), drop = FALSE]   # 与 pheno 一起按 GSM 排序，保证顺序一致

## ---------------------------------------------------------------- 3. 表型
pheno <- data.frame(
  gsm        = colnames(mat),
  sample_id  = colnames(mat),
  title      = title[match(colnames(mat), gsm)],
  donor_name = donor[match(colnames(mat), gsm)],
  region     = region[match(colnames(mat), gsm)],
  tbi_group  = tbi_group[match(colnames(mat), gsm)],
  dementia   = dem[match(colnames(mat), gsm)],
  sex        = sex[match(colnames(mat), gsm)],
  age_at_death = age[match(colnames(mat), gsm)],
  rin        = suppressWarnings(as.numeric(rin[match(colnames(mat), gsm)])),
  batch      = batch[match(colnames(mat), gsm)],
  braak      = braak[match(colnames(mat), gsm)],
  cerad      = cerad[match(colnames(mat), gsm)],
  apo_e4     = apo_e4[match(colnames(mat), gsm)],
  num_tbi    = ntbi[match(colnames(mat), gsm)],
  stringsAsFactors = FALSE
)
rownames(pheno) <- pheno$gsm

## ---------------------------------------------------------------- 4. 选分析单元
if (REGION != "ALL") {
  keep <- pheno$region == REGION
  if (sum(keep) == 0) stop("脑区 ", REGION, " 无样本。可选：",
                           paste(sort(unique(pheno$region)), collapse = "/"))
  mat   <- mat[, keep, drop = FALSE]
  pheno <- pheno[keep, , drop = FALSE]
}
message("[loader] 分析单元 region=", REGION, " → n=", ncol(mat),
        "（tbi ", sum(pheno$tbi_group == "tbi"), " / control ",
        sum(pheno$tbi_group == "control"), "），供者 ",
        length(unique(pheno$donor_name)))

## ---------------------------------------------------------------- 5. QC（§8.1）
qc <- list(
  gse = GSE, variant = VARIANT, region = REGION,
  n_samples = ncol(mat), n_genes = nrow(mat),
  n_tbi = sum(pheno$tbi_group == "tbi"),
  n_control = sum(pheno$tbi_group == "control"),
  n_donors = length(unique(pheno$donor_name)),
  samples_per_donor_max = max(table(pheno$donor_name)),
  n_batches = length(unique(pheno$batch)),
  min_batch_size = min(table(pheno$batch)),
  pure_batches = sum(tapply(pheno$tbi_group, pheno$batch,
                            function(z) length(unique(z)) == 1L)),
  fpkm_min = min(mat, na.rm = TRUE), fpkm_max = max(mat, na.rm = TRUE),
  fpkm_median = median(mat, na.rm = TRUE),
  n_na = sum(is.na(mat)),
  n_zero_genes = sum(rowSums(mat, na.rm = TRUE) == 0),
  scale_inferred = if (max(mat, na.rm = TRUE) > 50) "linear→log2(x+1) by 03" else "log_scale",
  colnames_equal_rownames_pheno = identical(colnames(mat), rownames(pheno))
)
for (k in names(qc)) message(sprintf("  [QC] %-32s %s", k, paste(qc[[k]], collapse = ",")))
if (!qc$colnames_equal_rownames_pheno) stop("矩阵列名与表型行名未对齐——03 会算错")
if (qc$n_na > 0) message("[loader] 警告：矩阵含 ", qc$n_na, " 个 NA")
if (qc$pure_batches > 0) message("[loader] 警告：存在 ", qc$pure_batches,
                                 " 个只含单一分组的批次（ComBat 会吞信号）")

## ---------------------------------------------------------------- 6. 落盘
ensure <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)
ensure(dirname(OUT))

eset <- Biobase::ExpressionSet(
  assayData = mat,
  phenoData = Biobase::AnnotatedDataFrame(pheno),
  annotation = NA_character_
)
saveRDS(list(eset), OUT)
message("[loader] 缓存已写入：", OUT, "（", file.size(OUT), " 字节）")

## QC 表与样本表一并留痕（审计用）
ensure("outputs")
write.csv(pheno, file.path("outputs", paste0("qc_samples_", GSE, "_", REGION, "_", VARIANT, ".csv")),
          row.names = FALSE)
write.csv(as.data.frame(qc, stringsAsFactors = FALSE),
          file.path("outputs", paste0("qc_summary_", GSE, "_", REGION, "_", VARIANT, ".csv")),
          row.names = FALSE)
message("[loader] 完成。manifest 需设：group_col=tbi_group, control_level=control, ",
        "case_level=tbi, patient_col=donor_name")
