## =====================================================================
## 08_determinants.R —— 可重复性的影响因素（对应手稿 §2.10）
##
## 层级得分是【有界比例】，所以用 beta 回归（logit 连接）；
## 零膨胀时用两部分模型（logistic 判有无 + beta 回归拟合正值部分）。
## 这条规则在预注册里就固定下来了，不是看到数据再决定的。
## 研究层面的聚类用聚类稳健标准误处理。
## 目标：给出【效应量的实际量级】，不是显著性检验。
## =====================================================================

## ---- 有界比例的常用变换（处理 0 与 1）----
## Smithson & Verkuilen 变换，n = 样本量；规则固定
transform_bounded <- function(y, eps_rule = "sv") {
  n <- length(y)
  (y * (n - 1) + 0.5) / n
}

## ---- 两部分模型的判据：零（或一）的比例超过阈值就用两部分 ----
needs_two_part <- function(y, cut = 0.10) {
  mean(y <= 0) > cut || mean(y >= 1) > cut
}

## ---- 聚类稳健标准误（能算就算，算不出就如实标注为朴素 SE）----
cluster_robust_se <- function(fit, cluster) {
  if (!requireNamespace("sandwich", quietly = TRUE) ||
      !requireNamespace("lmtest", quietly = TRUE)) {
    return(list(vcov = stats::vcov(fit), method = "naive_sandwich_missing"))
  }
  v <- try(sandwich::vcovCL(fit, cluster = cluster, type = "HC1"), silent = TRUE)
  if (inherits(v, "try-error")) {
    return(list(vcov = stats::vcov(fit), method = "naive_vcovcl_failed"))
  }
  list(vcov = v, method = "cluster_robust_HC1")
}

## ---- 单个结局变量的拟合 ----
fit_determinant <- function(y, X, cluster = NULL, label = "outcome") {
  y <- as.numeric(y)
  ok <- is.finite(y) & stats::complete.cases(X)
  y <- y[ok]; X <- X[ok, , drop = FALSE]
  if (!is.null(cluster)) cluster <- cluster[ok]
  if (nrow(X) < 5L || ncol(X) == 0L) {
    return(data.frame(outcome = label, term = NA_character_, estimate = NA_real_,
                      ci_lo = NA_real_, ci_hi = NA_real_, model = "insufficient_data",
                      se_method = NA_character_, stringsAsFactors = FALSE))
  }

  res <- list()
  if (needs_two_part(y)) {
    ## 第一部分：有无（logistic，报 OR）
    ybin <- as.integer(y > 0)
    fit1 <- try(suppressWarnings(stats::glm(ybin ~ ., data = data.frame(ybin, X),
                                            family = stats::binomial())), silent = TRUE)
    if (!inherits(fit1, "try-error")) {
      cr <- cluster_robust_se(fit1, cluster)
      ct <- if (requireNamespace("lmtest", quietly = TRUE)) {
        try(lmtest::coeftest(fit1, vcov = cr$vcov), silent = TRUE)
      } else NULL
      co <- if (!is.null(ct) && !inherits(ct, "try-error")) ct else
              summary(fit1)$coefficients
      for (i in seq_len(nrow(co))) {
        res[[length(res) + 1L]] <- data.frame(
          outcome = label, term = rownames(co)[i],
          estimate = exp(co[i, 1]),
          ci_lo = exp(co[i, 1] - 1.96 * co[i, 2]),
          ci_hi = exp(co[i, 1] + 1.96 * co[i, 2]),
          model = "two_part_logistic_any_vs_none",
          se_method = cr$method, stringsAsFactors = FALSE)
      }
    }
    ## 第二部分：正值部分的 beta 回归
    ypos <- y[y > 0]
    Xpos <- X[y > 0, , drop = FALSE]
    clpos <- if (!is.null(cluster)) cluster[y > 0] else NULL
    if (requireNamespace("betareg", quietly = TRUE) && length(ypos) > 5L &&
        length(unique(ypos)) > 2L) {
      ## 注意：响应变量必须是 data 里的列名 yv，否则 betareg 找不到它
      fit2 <- try(suppressWarnings(betareg::betareg(yv ~ .,
                                                    data = data.frame(yv = transform_bounded(ypos),
                                                                      Xpos, check.names = FALSE))),
                  silent = TRUE)
      if (!inherits(fit2, "try-error")) {
        co <- summary(fit2)$coefficients$mean
        for (i in seq_len(nrow(co))) {
          res[[length(res) + 1L]] <- data.frame(
            outcome = label, term = rownames(co)[i],
            estimate = exp(co[i, 1]),
            ci_lo = exp(co[i, 1] - 1.96 * co[i, 2]),
            ci_hi = exp(co[i, 1] + 1.96 * co[i, 2]),
            model = "two_part_beta_positive_part",
            se_method = "model_based", stringsAsFactors = FALSE)
        }
      }
    }
  } else {
    ## 单一 beta 回归
    if (!requireNamespace("betareg", quietly = TRUE)) {
      return(data.frame(outcome = label, term = NA_character_, estimate = NA_real_,
                        ci_lo = NA_real_, ci_hi = NA_real_,
                        model = "betareg_not_installed", se_method = NA_character_,
                        stringsAsFactors = FALSE))
    }
    yv <- transform_bounded(y)
    fit <- try(suppressWarnings(betareg::betareg(yv ~ .,
                                                 data = data.frame(yv, X))), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(data.frame(outcome = label, term = NA_character_, estimate = NA_real_,
                        ci_lo = NA_real_, ci_hi = NA_real_,
                        model = "betareg_failed", se_method = NA_character_,
                        stringsAsFactors = FALSE))
    }
    cr <- cluster_robust_se(fit, cluster)
    ct <- if (requireNamespace("lmtest", quietly = TRUE)) {
      try(lmtest::coeftest(fit, vcov = cr$vcov), silent = TRUE)
    } else NULL
    co <- if (!is.null(ct) && !inherits(ct, "try-error")) ct else
            summary(fit)$coefficients$mean
    for (i in seq_len(nrow(co))) {
      res[[length(res) + 1L]] <- data.frame(
        outcome = label, term = rownames(co)[i],
        estimate = exp(co[i, 1]),
        ci_lo = exp(co[i, 1] - 1.96 * co[i, 2]),
        ci_hi = exp(co[i, 1] + 1.96 * co[i, 2]),
        model = "beta_regression_logit_link",
        se_method = cr$method, stringsAsFactors = FALSE)
    }
  }
  if (length(res) == 0L) return(NULL)
  do.call(rbind, res)
}

## ---- 从清单构造预测变量矩阵（§2.10 列出的全部预测因子）----
build_predictors <- function(manifest) {
  m <- manifest
  getc <- function(nm) if (nm %in% names(m)) m[[nm]] else rep(NA, nrow(m))
  df <- data.frame(
    log_n_per_group  = log(as.numeric(getc("n_per_group"))),
    species          = factor(getc("species")),
    tissue           = factor(getc("tissue")),
    time_window      = factor(getc("time_window")),
    platform         = factor(getc("data_type")),
    origin_batch     = factor(getc("origin_batch_correction")),
    n_datasets_used  = as.numeric(getc("n_datasets_used")),
    ext_validation   = factor(getc("external_validation_reported")),
    hub_prespecified = factor(getc("hub_algorithm_prespecified")),
    pub_year         = as.numeric(getc("pub_year")),
    stringsAsFactors = FALSE
  )
  df
}

## ---- 全部结局 × 全部预测因子 ----
run_determinants <- function(outcomes, predictors, cluster = NULL) {
  rows <- lapply(names(outcomes), function(nm) {
    fit_determinant(outcomes[[nm]], predictors, cluster = cluster, label = nm)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) NULL else do.call(rbind, rows)
}
