source("00_config.R", encoding="UTF-8"); source("01_utils.R", encoding="UTF-8")
source("02_io.R", encoding="UTF-8"); source("03_preprocess.R", encoding="UTF-8")
source("04_de.R", encoding="UTF-8"); source("09_ortholog.R", encoding="UTF-8")
mf <- read.csv("data/datasets_template.csv", stringsAsFactors = FALSE)

build <- function(gse) {
  e <- readRDS(file.path("data", "raw", paste0(gse, ".rds")))[[1]]
  parts <- eset_to_parts(e)
  ph <- standardise_pheno(parts$pheno, "group", "control", "tbi",
                          patient_col = NULL, batch_col = NULL)
  pp <- preprocess_unified(parts$mat, ph, ann = NULL, dataset_id = gse)
  h  <- harmonize_gene_space(pp$mat, species = mf$species[match(gse, mf$gse)],
                             target = ORTHO_TARGET, dataset_id = gse)
  pp$mat <- h$mat
  pp$pheno <- pp$pheno[match(colnames(h$mat), pp$pheno$sample_id), , drop = FALSE]
  st <- run_de(pp, data_type = pp$data_type, paired = FALSE)
  cat(sprintf("%s: 样本 %d (case=%d/ctrl=%d)  映射后基因 %d  DEG=%d\n",
              gse, ncol(pp$mat), sum(pp$pheno$group == 1L), sum(pp$pheno$group == 0L),
              nrow(st), length(primary_de(st))))
  st
}

A <- build("GSE2871")
B <- build("GSE269748")
C <- build("GSE45997")

pair <- function(sA, sB, la, lb) {
  sa <- primary_de(sA); sb <- primary_de(sB)
  both <- intersect(sa, sb)
  if (length(both) == 0L) { cat(la, "vs", lb, "：无交集\n\n"); return(invisible(NULL)) }
  d <- data.frame(gene = both,
                  lfcA = sA$logFC[match(both, sA$gene)],
                  lfcB = sB$logFC[match(both, sB$gene)])
  d$same <- sign(d$lfcA) == sign(d$lfcB)
  cat("\n== ", la, " vs ", lb, "：双显著 n=", length(both),
      "  方向一致=", sum(d$same), " (", round(mean(d$same), 3), ")\n", sep = "")
  print(head(d[order(-abs(d$lfcA)), ], 15), row.names = FALSE)
  invisible(d)
}
d1 <- pair(A, B, "GSE2871", "GSE269748")
d2 <- pair(A, C, "GSE2871", "GSE45997")
d3 <- pair(C, B, "GSE45997", "GSE269748")

## ---- 标志基因方向体检（每个数据集）----
MK <- c("GFAP","LCN2","S100A8","S100A9","CCL2","TIMP1","CD68","AIF1","VIM","STAT3",
        "IL6","TNF","IL1B","C3","CXCL1","CXCL2","CXCL3","ARG1","MMP9","HMOX1",
        "PTGS2","SOCS3","CHI3L1","SERPINA3","OLFM1","NEFL","SNAP25","SYT1")
for (nm in c("GSE2871", "GSE45997", "GSE269748")) {
  s <- switch(nm, GSE2871 = A, GSE45997 = C, GSE269748 = B)
  ii <- match(MK, s$gene); ok <- !is.na(ii)
  cat("\n-- ", nm, " 标志基因 --\n", sep = "")
  print(data.frame(gene = MK[ok],
                   logFC = round(s$logFC[ii[ok]], 2),
                   fdr = signif(s$fdr[ii[ok]], 2)), row.names = FALSE)
}
