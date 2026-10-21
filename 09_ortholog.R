## =====================================================================
## 09_ortholog.R —— 跨物种基因标识统一（方案 A：Ensembl 一对一同源映射）
##
## 背景（为什么必须做）
##   纳入的 6 个数据集横跨 human / rat / mouse。同一基因在三个物种中的
##   symbol 大小写不同（人 TNF / 鼠 Tnf / 大鼠 Tnf），若直接按字符串求交集，
##   跨物种对的"重叠"会被人为压到 0，污染 L4/L5/L6/L7 的所有数字。
##
## 方案 A 的定义（2026-10-21 定案）
##   以 **人（human）** 为目标基因空间；非人数据集通过 Ensembl BioMart 的
##   ortholog_one2one 关系把基因 symbol 重命名为人的同源基因 symbol。
##   只保留严格一对一（两个方向都不允许重复），其余基因丢弃并记入审计表。
##
## 为什么绕开 biomaRt::getLDS（实测结论，勿改回）
##   1) getLDS 走双 mart，实测返回 HTTP 405 Method Not Allowed；
##   2) biomaRt 会把请求重定向到它 pin 住的归档站（本次为
##      jun2026.archive.ensembl.org），且 curl 超时只有 ~10 s，必超时；
##   3) 因此这里直接打 BioMart 的 martservice 端点（跟随 308 重定向、
##      超时 240 s），一次成功，并把结果缓存成 rds。
##
## 可复现性
##   映射表一旦抓取即缓存到 data/processed/ortholog_<from>_to_<to>.rds，
##   同时输出 outputs/ortholog_map_<from>_to_<to>.tsv 作为透明化附件。
##   缓存命中后不再联网 —— 后续任何一次重跑都用同一张表，结论不受
##   Ensembl 版本漂移影响。审计表记录 Ensembl release 与抓取时间。
## =====================================================================

ORTHO_TARGET       <- "human"      # 目标基因空间
BIOMART_URL        <- "https://www.ensembl.org/biomart/martservice"
BIOMART_TIMEOUT    <- 240L         # 秒；实测 ~16–25 s，留足余量
ORTHO_ONE2ONE      <- "ortholog_one2one"

## ---- 物种名 → Ensembl dataset 前缀 ----
.species_key <- function(s) {
  s <- tolower(trimws(as.character(s)[1]))
  switch(s,
    "human" = , "homo_sapiens" = , "hsapiens" = , "hs" = , "hsa" = "hsapiens",
    "mouse" = , "mus_musculus" = , "mmusculus" = , "mm"  = , "mmu" = "mmusculus",
    "rat"   = , "rattus_norvegicus" = , "rnorvegicus" = , "rn" = , "rno" = "rnorvegicus",
    stop("未知物种：", s, "（已支持：human / mouse / rat）"))
}

## ---- 组装 BioMart XML 查询：取 from→to 的同源基因名与同源类型 ----
.bm_ortholog_query <- function(from_key, to_key) {
  sprintf(paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
    '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" ',
    'count="" datasetConfigVersion="0.6">',
    '<Dataset name="%s_gene_ensembl" interface="default">',
    '<Attribute name="external_gene_name"/>',
    '<Attribute name="%s_homolog_associated_gene_name"/>',
    '<Attribute name="%s_homolog_orthology_type"/>',
    '</Dataset></Query>'), from_key, to_key, to_key)
}

## ---- 抓取（带缓存）一对一同源映射表 ----
## 返回 list(map = data.frame(from_gene, to_gene), release, url, fetched_at,
##           n_raw, n_one2one, from, to)
fetch_orthologs <- function(from = "mouse", to = ORTHO_TARGET,
                            cache = TRUE, cache_dir = "data/processed",
                            timeout = BIOMART_TIMEOUT) {
  fk <- .species_key(from); tk <- .species_key(to)

  ## 同源空间：无需映射
  if (identical(fk, tk)) {
    return(list(map = NULL, release = NA_character_, url = NA_character_,
                fetched_at = NA_character_, from = fk, to = tk,
                n_raw = 0L, n_one2one = 0L))
  }

  fpath <- file.path(cache_dir, sprintf("ortholog_%s_to_%s.rds", fk, tk))
  if (isTRUE(cache) && file.exists(fpath)) {
    obj <- readRDS(fpath)
    message(sprintf("[ortholog] 命中缓存 %s→%s：%d 对（Ensembl %s，抓取于 %s）",
                    fk, tk, nrow(obj$map), obj$release, obj$fetched_at))
    return(obj)
  }

  if (!requireNamespace("curl", quietly = TRUE)) {
    stop("需要 curl 包才能抓取 Ensembl 同源表：install.packages('curl')")
  }

  url <- paste0(BIOMART_URL, "?query=",
                utils::URLencode(.bm_ortholog_query(fk, tk), reserved = TRUE))
  message("[ortholog] 抓取 ", fk, "→", tk, "（最长 ", timeout, " s）…")
  h  <- curl::new_handle(followlocation = TRUE, timeout = timeout)
  rp <- curl::curl_fetch_memory(url, handle = h)
  if (rp$status_code != 200L) {
    stop(sprintf("BioMart HTTP %d（%s）", rp$status_code, substr(rawToChar(rp$content), 1, 200)))
  }
  body <- rawToChar(rp$content)
  eff  <- rp$url
  if (grepl("^\\s*(Query )?ERROR", body)) {
    stop("BioMart 返回错误：", substr(body, 1, 300))
  }

  ## 从最终 URL 反推 Ensembl release（如 jun2026.archive.ensembl.org → jun2026）
  rel <- regmatches(eff, regexpr("//[^./]+\\.archive\\.ensembl\\.org", eff))
  rel <- if (length(rel)) sub("//", "", sub("\\.archive\\.ensembl\\.org", "", rel)) else "www_latest"

  d <- utils::read.delim(text = body, header = FALSE, sep = "\t", quote = "",
                         comment.char = "", stringsAsFactors = FALSE,
                         col.names = c("from_gene", "to_gene", "orthology_type"))
  n_raw <- nrow(d)

  d <- d[!is.na(d$from_gene) & !is.na(d$to_gene) &
         d$from_gene != "" & d$to_gene != "", , drop = FALSE]
  d <- d[d$orthology_type == ORTHO_ONE2ONE, , drop = FALSE]

  ## 严格一对一：from 方向、to 方向都不允许重复
  d <- d[!(duplicated(d$from_gene) | duplicated(d$from_gene, fromLast = TRUE)), , drop = FALSE]
  d <- d[!(duplicated(d$to_gene)   | duplicated(d$to_gene,   fromLast = TRUE)), , drop = FALSE]

  map <- data.frame(from_gene = d$from_gene, to_gene = d$to_gene,
                    stringsAsFactors = FALSE)
  obj <- list(map = map, release = rel, url = eff,
              fetched_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              from = fk, to = tk, n_raw = n_raw, n_one2one = nrow(map))

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, fpath)
  message(sprintf("[ortholog] %s→%s：原始 %d 行 → 严格一对一 %d 对（Ensembl %s）",
                  fk, tk, n_raw, nrow(map), rel))
  obj
}

## ---- 把表达矩阵的行名映射到目标基因空间 ----
## 返回 list(mat, status, species, target, n_in, n_mapped,
##           n_unmapped_dropped, n_collapsed, ensembl_release, dataset_id)
## status: native_target_space / ok / failed_<原因>
## 失败时 mat = NULL —— 调用方必须把该数据集从跨数据集层级中剔除，
## 绝不允许"静默退回"原基因空间（那等于把污染的数字当结果用）。
harmonize_gene_space <- function(mat, species, target = ORTHO_TARGET,
                                 dataset_id = "dataset", cache = TRUE) {
  n_in <- nrow(mat)
  sp   <- tolower(trimws(as.character(species)[1]))
  base <- function(status, ...) list(
    mat = NULL, status = status, species = sp, target = target,
    n_in = n_in, n_mapped = 0L, n_unmapped_dropped = 0L, n_collapsed = 0L,
    ensembl_release = NA_character_, dataset_id = dataset_id, ...)

  if (is.na(sp) || sp == "") return(base("failed_unknown_species"))
  fk <- try(.species_key(sp), silent = TRUE)
  if (inherits(fk, "try-error")) return(base("failed_unknown_species"))

  ## 本来就人的基因空间：原样通过（但仍记一行审计，便于核对）
  if (identical(fk, .species_key(target))) {
    return(list(mat = mat, status = "native_target_space", species = sp,
                target = target, n_in = n_in, n_mapped = n_in,
                n_unmapped_dropped = 0L, n_collapsed = 0L,
                ensembl_release = NA_character_, dataset_id = dataset_id))
  }

  o <- try(fetch_orthologs(sp, target, cache = cache), silent = TRUE)
  if (inherits(o, "try-error")) {
    message("[ortholog][失败] ", dataset_id, "：", conditionMessage(attr(o, "condition")))
    return(base("failed_fetch", error = substr(conditionMessage(attr(o, "condition")), 1, 200)))
  }

  idx  <- match(rownames(mat), o$map$from_gene)
  keep <- which(!is.na(idx))
  if (length(keep) == 0L) {
    message("[ortholog][失败] ", dataset_id, "：没有任何基因能映射到 ", target, " 空间")
    return(base("failed_empty_map", ensembl_release = o$release))
  }

  m2   <- mat[keep, , drop = FALSE]
  newn <- o$map$to_gene[idx[keep]]
  n_unmapped <- n_in - length(keep)

  ## 多个源基因映射到同一目标 symbol 时，保留平均表达最高的那一行
  dupn <- unique(newn[duplicated(newn)])
  if (length(dupn)) {
    rm_ <- rowMeans(m2, na.rm = TRUE)
    pick <- unlist(lapply(dupn, function(s) {
      ii <- which(newn == s); ii[which.max(rm_[ii])]
    }))
    keep2 <- sort(c(which(!newn %in% dupn), as.integer(pick)))
    m2 <- m2[keep2, , drop = FALSE]; newn <- newn[keep2]
  }
  rownames(m2) <- newn

  message(sprintf("[ortholog] %s（%s→%s）：%d 行 → %d 行（未映射丢弃 %d，重复折叠 %d）",
                  dataset_id, sp, target, n_in, nrow(m2), n_unmapped, length(dupn)))

  list(mat = m2, status = "ok", species = sp, target = target,
       n_in = n_in, n_mapped = nrow(m2), n_unmapped_dropped = n_unmapped,
       n_collapsed = length(dupn), ensembl_release = o$release,
       dataset_id = dataset_id)
}

## ---- 审计行：一条数据集一行，必须随结果一起交付 ----
ortholog_audit_row <- function(h) {
  data.frame(
    dataset_id          = h$dataset_id,
    species             = h$species,
    target_space        = h$target,
    status              = h$status,
    n_genes_in          = h$n_in,
    n_genes_mapped      = h$n_mapped,
    n_unmapped_dropped  = h$n_unmapped_dropped,
    n_collapsed         = h$n_collapsed,
    ensembl_release     = if (is.null(h$ensembl_release)) NA_character_ else h$ensembl_release,
    stringsAsFactors = FALSE)
}
