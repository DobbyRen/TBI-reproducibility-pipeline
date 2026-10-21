## =====================================================================
## make_gmt.R —— 生成 L6 需要的 MSigDB 基因集 .gmt 文件
##
## 实测结论（2026-10-21，勿改回旧路径）
##   1) MSigDB 官网 gsea-msigdb.org 下载【需要登录】，未登录时返回的其实是
##      一个 4425 字节的 HTML 登录页（HTTP 200，很容易被误判为成功）。
##   2) msigdbr 包的数据托管在 zenodo，实测下载速率约 10 KB/s，
##      100+ MB 的压缩包基本下不动。
##   3) Broad 的开放镜像 data.broadinstitute.org 可直接取 v7.5.1 的
##      symbols.gmt（28 MB，实测 5.7 MB/s），无需登录，且一次下载就同时
##      包含 Hallmark / KEGG / GO_BP —— 本脚本采此路径。
##
## 用法：Rscript make_gmt.R
## 产出：data/gmt/msigdb_{Hallmark,KEGG,GOBP}_v7.5.1.gmt + README_provenance.md
## =====================================================================

MSIGDB_URL  <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/msigdb.v7.5.1.symbols.gmt"
MSIGDB_VER  <- "7.5.1"
OUT_DIR     <- "data/gmt"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- 1. 下载（本地已有同名缓存则跳过）----
cache <- file.path(tempdir(), paste0("msigdb.v", MSIGDB_VER, ".symbols.gmt"))
if (!file.exists(cache)) {
  cat("[gmt] 下载 MSigDB v", MSIGDB_VER, " …\n", sep = "")
  if (!requireNamespace("curl", quietly = TRUE)) stop("需要 curl 包")
  h <- curl::new_handle(followlocation = TRUE, timeout = 600L)
  r <- curl::curl_fetch_disk(MSIGDB_URL, cache, handle = h)
  if (r$status_code != 200L) stop("下载失败，HTTP ", r$status_code)
} else {
  cat("[gmt] 使用缓存：", cache, "\n")
}
cat("[gmt] 文件大小：", round(file.size(cache) / 1024 / 1024, 2), " MB\n")

## 校验：登录页是 HTML，真 gmt 是 ASCII 文本且首字段不含 '<'
head1 <- readLines(cache, n = 1, warn = FALSE)
if (grepl("<", head1, fixed = TRUE)) stop("下载到的是 HTML 页面而非 .gmt，请检查是否需要登录")

## ---- 2. 按集合前缀拆分 ----
ln <- readLines(cache, warn = FALSE)
pick <- function(prefix) ln[grepl(paste0("^", prefix), ln)]
write_part <- function(lines, path) {
  writeLines(lines, path)
  length(lines)
}

parts <- c(Hallmark = "HALLMARK_", KEGG = "KEGG_", GOBP = "GOBP_")
n <- integer(0)
for (nm in names(parts)) {
  out <- file.path(OUT_DIR, sprintf("msigdb_%s_v%s.gmt", nm, MSIGDB_VER))
  n[nm] <- write_part(pick(parts[[nm]]), out)
  cat(sprintf("[gmt] %-9s %5d 个集合 → %s\n", nm, n[nm], out))
}

## ---- 3. 来源台账（投稿时要能说清版本与获取途径）----
prov <- c(
  "# 基因集来源（data/gmt/）",
  "",
  sprintf("- 数据库：MSigDB v%s（Human, symbols）", MSIGDB_VER),
  sprintf("- 下载地址：%s", MSIGDB_URL),
  sprintf("- 本地文件 md5：%s", tools::md5sum(cache)),
  sprintf("- 拆分时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("- Hallmark %d 个集合 / KEGG %d 个 / GO_BP %d 个", n["Hallmark"], n["KEGG"], n["GOBP"]),
  "",
  "说明：MSigDB 官网下载需登录；Broad 开放镜像无需登录且为同一份数据。",
  "所有非人数据集在 L6 之前已通过 09_ortholog.R 映射到人类基因空间，",
  "因此只用人类基因集即可覆盖全部数据集。"
)
writeLines(prov, file.path(OUT_DIR, "README_provenance.md"))

cat("[gmt] 完成。下一步：Rscript run_all.R synthetic（L6 将不再被跳过）\n")
