## =====================================================================
## 01_utils.R —— 通用工具函数（无外部依赖，纯 base R）
## =====================================================================

## ---- 目录 ----
ensure_dirs <- function(...) {
  for (d in c(...)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
  invisible(NULL)
}

## ---- 集合度量：L2 用 ----
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  if (union == 0) return(NA_real_)
  inter / union
}

overlap_coef <- function(a, b) {
  a <- unique(a); b <- unique(b)
  inter <- length(intersect(a, b))
  denom <- min(length(a), length(b))
  if (denom == 0) return(NA_real_)
  inter / denom
}

## 超几何检验：观测重叠是否超过随机期望（背景 = universe）
hyper_p <- function(a, b, universe, alternative = "greater") {
  a <- unique(a); b <- unique(b)
  n_uni <- length(unique(universe))
  k <- length(a)                      # A 的规模
  m <- length(b)                      # 背景中"命中"的规模（B）
  q <- length(intersect(a, b))
  if (n_uni == 0 || k == 0 || m == 0) return(NA_real_)
  phyper(q - 1, m, n_uni - m, k, lower.tail = FALSE, log.p = FALSE)
}

## ---- 自助法置信区间 ----
boot_ci_prop <- function(x, B = FROZEN$boot_b, conf = FROZEN$boot_ci, seed = NULL) {
  x <- as.integer(as.logical(x))
  n <- length(x)
  if (n == 0) return(c(est = NA_real_, lo = NA_real_, hi = NA_real_, n = 0L))
  if (!is.null(seed)) set.seed(seed)
  est <- mean(x)
  bs <- replicate(B, mean(sample(x, n, replace = TRUE)))
  q <- as.numeric(quantile(bs, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE))
  c(est = est, lo = q[1], hi = q[2], n = n)
}

boot_ci_diff <- function(x, y, B = FROZEN$boot_b, conf = FROZEN$boot_ci, seed = NULL) {
  ## x、y 为 0/1 向量，返回 rate_x - rate_y 的点估计与自助 CI
  if (!is.null(seed)) set.seed(seed)
  nx <- length(x); ny <- length(y)
  if (nx == 0 || ny == 0) return(c(diff = NA_real_, lo = NA_real_, hi = NA_real_))
  d_hat <- mean(x) - mean(y)
  bs <- replicate(B, mean(sample(as.integer(x), nx, TRUE)) -
                     mean(sample(as.integer(y), ny, TRUE)))
  q <- as.numeric(quantile(bs, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE))
  c(diff = d_hat, lo = q[1], hi = q[2])
}

## ---- top-k 交集：L3 用 ----
topk_overlap <- function(rankA, rankB, k) {
  a <- head(unique(rankA), k)
  b <- head(unique(rankB), k)
  length(intersect(a, b)) / k
}

## ---- 排名：FDR 升序，FDR 相同时按 |log2FC| 降序（预注册锁定）----
rank_genes <- function(stats) {
  stats <- stats[order(stats$fdr, -abs(stats$logFC), stats$gene), ]
  stats$gene
}

## ---- 主阈值筛选 ----
primary_de <- function(stats, fdr = FROZEN$de_fdr_primary, lfc = FROZEN$de_lfc_primary) {
  stats$gene[stats$fdr < fdr & abs(stats$logFC) >= lfc]
}

## ---- 二分类性能度量（不依赖 pROC / PRROC）----
auc_mw <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]; y <- as.integer(y[ok])
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

pr_auc <- function(score, y) {
  ok <- is.finite(score) & !is.na(y)
  score <- score[ok]; y <- as.integer(y[ok])
  if (length(unique(y)) < 2) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  recall <- tp / max(sum(y == 1), 1)
  precis <- tp / pmax(tp + fp, 1)
  ## 梯形法求平均精度（baseline = 正例比例，作为参考一起返回）
  x <- c(0, recall); z <- c(1, precis)
  ap <- sum(diff(x) * (head(z, -1) + tail(z, -1)) / 2)
  c(ap = ap, baseline = mean(y))
}

brier <- function(p, y) {
  ok <- is.finite(p) & !is.na(y)
  mean((p[ok] - as.integer(y[ok])) ^ 2)
}

logloss <- function(p, y, eps = 1e-12) {
  ok <- is.finite(p) & !is.na(y)
  p <- pmin(pmax(p[ok], eps), 1 - eps)
  yy <- as.integer(y[ok])
  -mean(yy * log(p) + (1 - yy) * log(1 - p))
}

clip01 <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)

## ---- 校准：slope / intercept / 曲线 ----
calibration_fit <- function(p, y) {
  ok <- is.finite(p) & !is.na(y) & length(unique(y)) > 1
  p <- clip01(p[ok]); y <- as.integer(y[ok])
  lp <- log(p / (1 - p))
  fit <- try(suppressWarnings(glm(y ~ lp, family = binomial())), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(slope = NA_real_, intercept = NA_real_, curve = NULL))
  }
  co <- coef(fit)
  ## 分 10 组绘制校准曲线用（论文 Fig. 6）
  grp <- cut(p, breaks = quantile(p, probs = seq(0, 1, 0.1), na.rm = TRUE),
             include.lowest = TRUE, labels = FALSE)
  agg <- data.frame(bin = grp, p = p, y = y)
  curve <- aggregate(cbind(p, y) ~ bin, data = agg, FUN = mean)
  list(slope = unname(co[2]), intercept = unname(co[1]), curve = curve)
}

## ---- 交叉验证折叠（按患者水平）----
make_cv_folds <- function(patient, y, folds = FROZEN$cv_folds, seed) {
  set.seed(seed)
  pats <- unique(patient)
  ## 按患者分层，保证每折两组的比例接近
  pat_y <- tapply(y, patient, function(v) as.integer(mean(v) > 0.5))
  idx <- split(sample(pats),
               cut(seq_along(pats), breaks = folds, labels = FALSE))
  lapply(idx, function(pp) which(patient %in% pp))
}

## ---- 输出辅助 ----
write_out <- function(x, name, dir = "outputs") {
  ensure_dirs(dir)
  f <- file.path(dir, paste0(name, ".csv"))
  utils::write.csv(x, f, row.names = FALSE)
  message("[out] ", f)
  invisible(f)
}
