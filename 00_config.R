## =====================================================================
## 00_config.R —— 冻结参数（与预注册方案逐条对应）
## 文章 I：TBI 转录组生物标志物声明的可重复性 —— 七级定量再分析
##
## 纪律：本文件在 OSF 预注册提交后不得修改。任何修改都必须通过
##       OSF 注册的「Update」记录，并同步写入协议 Appendix A
##       修订台账（Amendment log），两者逐条对应。
##       log_frozen() 会把本文件的 MD5 与全部参数写进运行日志，
##       用于事后证明"跑结果时用的就是这一套参数"。
## =====================================================================

FROZEN <- list(

  ## ---- 随机种子：全流程统一派生，禁止在函数内部另设 seed ----
  seed_master = 20260501L,

  ## ---- 差异表达：主分析阈值（预注册锁定，不可调）----
  de_fdr_primary = 0.05,
  de_lfc_primary = 0.50,

  ## ---- 差异表达：稳健性阈值网格 ----
  ## 注意：预注册正文写的是 |log2FC| ∈ {0.25, 1.00}；
  ##       手稿 §2.8 写的是 |log2FC| ∈ {0, 0.5, 1}。
  ##       此处取并集 {0, 0.25, 0.50, 1.00}。
  ##       主分析阈值（FDR<0.05 且 |log2FC|≥0.5）两者一致，故结论判据未变，
  ##       但这属于"扩充网格"，按 §9 惯例应记入修订台账 v1.1。
  de_fdr_grid = c(0.01, 0.05, 0.10),
  de_lfc_grid = c(0, 0.25, 0.50, 1.00),

  ## ---- 过滤：至少在一个组中 ≥20% 的样本表达量高于固定下限 ----
  filter_min_frac = 0.20,
  filter_floor    = list(log_scale = 4,      # log2 尺度：4 ≈ 线性 16
                         linear    = 16,
                         count     = 16),

  ## ---- 探针聚合 ----
  probe_agg_primary     = "maxmean",  # 同基因多探针：保留平均表达最高的探针
  probe_agg_sensitivity = "mean",     # 敏感性分析：取均值
  probe_drop_multi      = TRUE,       # 映射到多个基因的探针直接丢弃

  ## ---- 批次：仅当数据集显式提供批次变量时才做 ComBat，且保护分组变量 ----
  batch_method = "combat",

  ## ---- 内部双层拆分（REDE-2Fold 的 TBI 适配）----
  split_level    = "patient",  # 患者/动物水平，绝不能是样本水平
  split_prop     = 0.5,
  split_stratify = TRUE,

  ## ---- 负对照：规模匹配的随机抽样 ----
  negctrl_draws = 1000L,

  ## ---- L3 高排名重叠 ----
  topk_grid = c(25L, 50L, 100L),

  ## ---- L5 预先指定基因的确认 ----
  l5_lfc_large = 1.00,   # 大效应子集判据

  ## ---- L6 通路层面 ----
  enrich_fdr      = 0.05,
  enrich_dbs      = c("Hallmark", "KEGG", "GO_BP"),
  enrich_min_set  = 10L,
  enrich_max_set  = 500L,
  enrich_engine   = "auto",   # "clusterProfiler" | "ora_hyper" | "auto"

  ## ---- L7 建模与操作点 ----
  model_primary     = "glmnet_l2",
  model_sensitivity = c("rf", "gbm"),
  cv_folds          = 5L,
  cv_repeats        = 5L,
  op_rule_primary     = "youden",            # 主分析：Youden 最优点
  op_rule_sensitivity = "fixed_sens90",      # 敏感性：固定敏感度 90%
  op_target_sens      = 0.90,
  op_failure_cut      = 0.10,                # 锁定操作点下 sens 或 spec < 0.10 判为失效

  ## ---- 自助法与检验 ----
  boot_b        = 2000L,
  boot_ci       = 0.95,
  sign_test_alt = "greater",   # L4 符号检验：H0 方向一致概率 = 0.5

  ## ---- L4/L5 的方向要求 ----
  require_direction = TRUE     # 确认必须"显著 + 方向一致"，显著但反向不算
)

## ---------------------------------------------------------------------
## 参数指纹：把本文件的 MD5 与全部参数写进日志
## ---------------------------------------------------------------------
log_frozen <- function(config_path = "00_config.R", logdir = "outputs") {
  dir.create(logdir, showWarnings = FALSE, recursive = TRUE)
  md5 <- if (file.exists(config_path)) unname(tools::md5sum(config_path)) else NA_character_
  flat <- vapply(FROZEN,
                 function(x) paste(format(x, trim = TRUE), collapse = ", "),
                 character(1))
  lines <- c(
    "## FROZEN PARAMETERS — do not edit after pre-registration",
    paste0("timestamp   : ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("config_file : ", config_path),
    paste0("config_md5  : ", md5),
    paste0("R_version   : ", R.version.string),
    "",
    paste0("  ", names(FROZEN), " = ", flat)
  )
  writeLines(lines, file.path(logdir, "frozen_params.log"))
  message("[config] frozen parameters logged -> ", file.path(logdir, "frozen_params.log"))
  invisible(lines)
}

## ---------------------------------------------------------------------
## 派生种子：把字符标签稳定映射为整数种子，避免所有环节共用一个种子
## ---------------------------------------------------------------------
seed_for <- function(tag) {
  ii <- utf8ToInt(paste0("TBI::", tag))
  h <- sum(ii * seq_along(ii))
  as.integer((FROZEN$seed_master + h) %% 2147483000L)
}

## ---------------------------------------------------------------------
## 断言：任何试图在主分析之外改动主阈值的地方都应在入口调用本函数
## ---------------------------------------------------------------------
assert_primary_threshold <- function(fdr = FROZEN$de_fdr_primary,
                                     lfc = FROZEN$de_lfc_primary) {
  if (!isTRUE(all.equal(fdr, FROZEN$de_fdr_primary)) ||
      !isTRUE(all.equal(lfc, FROZEN$de_lfc_primary))) {
    warning("[config] 非主分析阈值（FDR=", fdr, ", |log2FC|=", lfc,
            "）——请确保这属于敏感性分析并已在论文中标注。")
  }
  invisible(TRUE)
}
