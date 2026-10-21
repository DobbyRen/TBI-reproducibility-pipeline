## =====================================================================
## 07_model_transfer.R —— L7 判别性能与锁定操作点迁移（对应手稿 §2.9）
##
## 全文最容易被审稿人抓的地方，三条铁律：
##  1. 在接触外部数据【之前】锁定：特征集、预处理、模型、超参数、决策阈值。
##     本文件的函数签名就是纪律——lock_operating_point() 只接受训练数据，
##     物理上拿不到外部数据的任何信息。
##  2. 阈值用训练集交叉验证的 out-of-fold 预测值确定，避免乐观偏倚。
##  3. 必须报告"高 AUC 但锁定操作点下敏感度或特异度接近零"的模型比例。
##     这是 H7 的靶心，也是常规报道里最容易看不见的失效模式。
## =====================================================================

## ---- 建模：主分析 logistic regression + L2（glmnet）----
build_glmnet <- function(x, y, patient = NULL, seed) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("需要 glmnet：install.packages('glmnet')")
  }
  x <- as.matrix(x); y <- as.integer(y)
  set.seed(seed)
  lambdas <- numeric(FROZEN$cv_repeats)
  cvms <- numeric(FROZEN$cv_repeats)
  best <- NULL
  for (r in seq_len(FROZEN$cv_repeats)) {
    foldid <- if (is.null(patient)) {
      sample(rep_len(seq_len(FROZEN$cv_folds), length(y)))
    } else {
      fo <- make_cv_folds(patient, y, folds = FROZEN$cv_folds,
                          seed = seed_for(paste0("cv", r)))
      fid <- integer(length(y))
      for (i in seq_along(fo)) fid[fo[[i]]] <- i
      fid
    }
    cv <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 0,
                            foldid = foldid, type.measure = "auc")
    lambdas[r] <- cv$lambda.min
    cvms[r] <- min(cv$cvm)
    if (is.null(best) || min(cv$cvm) < min(cvms, na.rm = TRUE)) best <- cv
  }
  lam <- stats::median(lambdas)
  fit <- glmnet::glmnet(x, y, family = "binomial", alpha = 0, lambda = lam)
  list(fit = fit, lambda = lam, cv_auc = mean(cvms, na.rm = TRUE), method = "glmnet_l2")
}

## ---- 敏感性分析：随机森林 ----
build_rf <- function(x, y, seed) {
  if (requireNamespace("ranger", quietly = TRUE)) {
    set.seed(seed)
    fit <- ranger::ranger(y = factor(y), x = as.data.frame(x),
                          probability = TRUE, importance = "none")
    return(list(fit = fit, method = "rf_ranger"))
  }
  if (requireNamespace("randomForest", quietly = TRUE)) {
    set.seed(seed)
    fit <- randomForest::randomForest(x = as.data.frame(x), y = factor(y))
    return(list(fit = fit, method = "rf_randomForest"))
  }
  message("[model] ranger / randomForest 均未安装，跳过随机森林")
  NULL
}

## ---- 敏感性分析：梯度提升 ----
build_gbm <- function(x, y, seed) {
  if (requireNamespace("xgboost", quietly = TRUE)) {
    set.seed(seed)
    fit <- xgboost::xgboost(data = as.matrix(x), label = as.integer(y),
                            nrounds = 100L, max_depth = 3L, eta = 0.1,
                            objective = "binary:logistic", verbose = 0)
    return(list(fit = fit, method = "gbm_xgboost"))
  }
  if (requireNamespace("gbm", quietly = TRUE)) {
    set.seed(seed)
    fit <- gbm::gbm(y ~ ., data = data.frame(y = as.integer(y), as.data.frame(x)),
                    distribution = "bernoulli", n.trees = 200L,
                    interaction.depth = 3L, shrinkage = 0.1, verbose = FALSE)
    return(list(fit = fit, method = "gbm_gbm"))
  }
  message("[model] xgboost / gbm 均未安装，跳过梯度提升")
  NULL
}

build_model <- function(x, y, patient = NULL,
                        method = FROZEN$model_primary, seed) {
  switch(method,
         glmnet_l2 = build_glmnet(x, y, patient, seed),
         rf        = build_rf(x, y, seed),
         gbm       = build_gbm(x, y, seed),
         stop("未知模型：", method))
}

## ---- 预测概率：统一接口 ----
predict_p <- function(model, x) {
  m <- model$method
  if (m == "glmnet_l2") {
    as.numeric(predict(model$fit, newx = as.matrix(x), s = model$lambda,
                       type = "response"))
  } else if (m == "rf_ranger") {
    predict(model$fit, data = as.data.frame(x))$predictions[, "1"]
  } else if (m == "rf_randomForest") {
    predict(model$fit, newdata = as.data.frame(x), type = "prob")[, "1"]
  } else if (m == "gbm_xgboost") {
    as.numeric(predict(model$fit, newdata = as.matrix(x)))
  } else if (m == "gbm_gbm") {
    plogis(predict(model$fit, newdata = as.data.frame(x), n.trees = 200L,
                   type = "link"))
  } else stop("未实现的预测接口：", m)
}

## ---- out-of-fold 预测（用于锁定操作点，避免用训练集内预测造成乐观）----
oof_predictions <- function(x, y, patient = NULL, seed) {
  fo <- if (is.null(patient)) {
    set.seed(seed); split(sample(seq_along(y)),
                          cut(seq_along(y), breaks = FROZEN$cv_folds, labels = FALSE))
  } else {
    make_cv_folds(patient, y, folds = FROZEN$cv_folds, seed = seed)
  }
  p <- numeric(length(y))
  for (i in seq_along(fo)) {
    tr <- setdiff(seq_along(y), fo[[i]])
    md <- build_glmnet(x[tr, , drop = FALSE], y[tr], patient = patient[tr],
                       seed = seed_for(paste0("oof", i)))
    p[fo[[i]]] <- predict_p(md, x[fo[[i]], , drop = FALSE])
  }
  p
}

## ---- 锁定操作点：只吃训练数据 ----
## rule = "youden"         : 最大化 sens + spec - 1
## rule = "fixed_sens90"   : 在满足敏感度 ≥90% 的阈值中取最大（特异度最高）
lock_operating_point <- function(score_train, y_train,
                                 rule = FROZEN$op_rule_primary) {
  ok <- is.finite(score_train) & !is.na(y_train)
  s <- score_train[ok]; y <- as.integer(y_train[ok])
  stopifnot(all(y %in% c(0L, 1L)))
  cuts <- sort(unique(s))
  if (length(cuts) == 0L) stop("训练集预测值全为 NA，无法锁定操作点")
  perf <- vapply(cuts, function(t) {
    pred <- as.integer(s >= t)
    tp <- sum(pred == 1L & y == 1L); fn <- sum(pred == 0L & y == 1L)
    tn <- sum(pred == 0L & y == 0L); fp <- sum(pred == 1L & y == 0L)
    c(sens = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
      spec = if (tn + fp > 0) tn / (tn + fp) else NA_real_)
  }, numeric(2))
  if (rule == "youden") {
    jd <- perf["sens", ] + perf["spec", ] - 1
    thr <- cuts[which.max(jd)]
  } else if (rule == "fixed_sens90") {
    ok2 <- which(perf["sens", ] >= FROZEN$op_target_sens)
    thr <- if (length(ok2)) max(cuts[ok2]) else min(cuts)
  } else stop("未知操作点规则：", rule)
  idx <- which(cuts == thr)[1]
  list(
    threshold = thr, rule = rule,
    train_sens = perf["sens", idx], train_spec = perf["spec", idx],
    locked_at = Sys.time(),
    note = "阈值由训练集 out-of-fold 预测确定，此后不再调整"
  )
}

## ---- 外部评估：锁定操作点下的敏感度/特异度/预测值 ----
evaluate_locked <- function(threshold, score_ext, y_ext, prevalence = NULL) {
  ok <- is.finite(score_ext) & !is.na(y_ext)
  s <- score_ext[ok]; y <- as.integer(y_ext[ok])
  pred <- as.integer(s >= threshold)
  tp <- sum(pred == 1L & y == 1L); fn <- sum(pred == 0L & y == 1L)
  tn <- sum(pred == 0L & y == 0L); fp <- sum(pred == 1L & y == 0L)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  obs_prev <- mean(y)
  ppv_obs <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  npv_obs <- if (tn + fn > 0) tn / (tn + fn) else NA_real_

  ## 按目标患病率换算（论文需给出 PPV/NPV 而非只给 sens/spec）
  ppv_tgt <- npv_tgt <- NA_real_
  if (!is.null(prevalence) && prevalence > 0 && prevalence < 1) {
    num <- sens * prevalence
    ppv_tgt <- num / (num + (1 - spec) * (1 - prevalence))
    num2 <- spec * (1 - prevalence)
    npv_tgt <- num2 / (num2 + (1 - sens) * prevalence)
  }

  ## 判别与校准
  p_auc <- auc_mw(s, y)
  pr    <- pr_auc(s, y)
  probs <- clip01(s)
  cal   <- calibration_fit(probs, y)

  data.frame(
    n_ext = length(y), observed_prevalence = obs_prev,
    auc = p_auc, pr_auc = pr[["ap"]], pr_baseline = pr[["baseline"]],
    brier = brier(probs, y), log_loss = logloss(probs, y),
    calib_slope = cal$slope, calib_intercept = cal$intercept,
    locked_threshold = threshold,
    sens = sens, spec = spec,
    ppv_observed = ppv_obs, npv_observed = npv_obs,
    prevalence_target = if (is.null(prevalence)) NA_real_ else prevalence,
    ppv_target = ppv_tgt, npv_target = npv_tgt,
    stringsAsFactors = FALSE
  )
}

## ---- 跨物种迁移：仅保留一对一同源基因（可选）----
ortholog_map <- function(genes, from = "mmusculus", to = "hsapiens") {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    message("[ortholog] biomaRt 未安装，跳过同源映射")
    return(list(genes = genes, n_discarded = 0L))
  }
  m <- try(biomaRt::getLDS(attributes = "external_gene_name",
                           filters = "external_gene_name", values = genes,
                           mart = biomaRt::useMart("ensembl", dataset = paste0(from, "_gene_ensembl")),
                           attributesL = "external_gene_name",
                           martL = biomaRt::useMart("ensembl", dataset = paste0(to, "_gene_ensembl"))),
           silent = TRUE)
  if (inherits(m, "try-error")) {
    message("[ortholog] 同源映射失败，退回原基因集")
    return(list(genes = genes, n_discarded = 0L))
  }
  mapped <- unique(m[[2]])
  list(genes = intersect(genes, mapped), n_discarded = length(genes) - length(intersect(genes, mapped)))
}

## ---- 一次完整的 L7 审计：训练 → 锁定 → 外部评估 ----
run_l7 <- function(x_train, y_train, patient_train,
                   x_ext, y_ext,
                   method = FROZEN$model_primary,
                   rule = FROZEN$op_rule_primary,
                   prevalence = NULL, seed, label = "model") {
  ## 1) 在训练集内做 CV，得到 out-of-fold 预测
  p_oof <- oof_predictions(x_train, y_train, patient_train, seed = seed)

  ## 2) 在接触外部数据之前锁定阈值
  op <- lock_operating_point(p_oof, y_train, rule = rule)

  ## 3) 用全部训练数据拟合最终模型
  md <- build_model(x_train, y_train, patient_train, method = method, seed = seed)

  ## 4) 外部评估（此后再无任何调参）
  p_ext <- predict_p(md, x_ext)
  ev <- evaluate_locked(op$threshold, p_ext, y_ext, prevalence = prevalence)

  ev$model_id <- label
  ev$method <- method
  ev$op_rule <- rule
  ev$train_cv_auc <- auc_mw(p_oof, y_train)
  ev$train_sens_at_lock <- op$train_sens
  ev$train_spec_at_lock <- op$train_spec
  ev$auc_drop <- ev$train_cv_auc - ev$auc
  ev$locked_point_failed <- (is.na(ev$sens) | ev$sens < FROZEN$op_failure_cut) |
                            (is.na(ev$spec) | ev$spec < FROZEN$op_failure_cut)
  ev
}

## ---- 失效比例：H7 的核心输出 ----
summarise_l7 <- function(tab) {
  n <- nrow(tab)
  data.frame(
    n_models = n,
    median_auc = stats::median(tab$auc, na.rm = TRUE),
    median_auc_drop = stats::median(tab$auc_drop, na.rm = TRUE),
    prop_locked_failed = if (n > 0) mean(tab$locked_point_failed) else NA_real_,
    prop_sens_lt_cut = if (n > 0) mean(tab$sens < FROZEN$op_failure_cut, na.rm = TRUE) else NA_real_,
    prop_spec_lt_cut = if (n > 0) mean(tab$spec < FROZEN$op_failure_cut, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}
