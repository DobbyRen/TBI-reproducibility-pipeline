#!/usr/bin/env Rscript
# 生成文章 I 的 7 张结果图（PDF 矢量，经 sips 转 PNG）
# 数据源全部来自 TBI_repro_pipeline/outputs/*.csv 与仲裁后合并表统计
# 只画图，不做任何统计推断；所有数字均读自既有结果文件。

OUT    <- "/Users/renqingqing/WorkBuddy/2026-09-30-10-42-09/outputs/manuscript/figures"
OUTPDF <- file.path(OUT, "pdf")
dir.create(OUTPDF, recursive = TRUE, showWarnings = FALSE)
setwd("/Users/renqingqing/WorkBuddy/2026-09-30-10-42-09/TBI_repro_pipeline")

rd <- function(f) read.csv(file.path("outputs", f), stringsAsFactors = FALSE)

de   <- rd("L1_de_burden.csv")
l2   <- rd("L2_membership.csv")
l3   <- rd("L3_topk_overlap.csv")
l4   <- rd("L4_concordance.csv")
stab <- rd("internal_stability.csv")
l5d  <- rd("L5_literature_by_dataset.csv")
l5s  <- rd("L5_literature_by_study.csv")

COL <- c("GSE104687"="#2C5F8A","GSE160763"="#8E44AD","GSE269748"="#B03A2E",
         "GSE2871"="#1B7F5C","GSE45997"="#D68910","GSE89866"="#566573")
col_main <- "#2C5F8A"; col_acc <- "#B03A2E"; col_grey <- "#8C8C8C"
col_ok <- "#1B7F5C"

pngf <- function(name, w, h) {
  # 本机 png/cairo/quartz 设备不可用，用 pdf() 出矢量图，sips 转 PNG
  pdf(file.path(OUTPDF, paste0(name, ".pdf")), width = w, height = h,
      bg = "white", pointsize = 9, useDingbats = FALSE)
}

short <- function(x) sub("^GSE", "", x)
wil <- function(x, n) {
  z <- 1.959964; p <- x / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  cbind(pmax(0, ctr - hw), pmin(1, ctr + hw))
}

## ---------------- Fig 1. PRISMA 2020 flow diagram ----------------
pngf("Fig1", 8, 10)
op <- par(mar = rep(0.4, 4)); on.exit(par(op))
plot.new(); plot.window(c(0, 100), c(0, 128))
boxx <- function(x, y, w, h, lab, sub = NULL, lwd = 1, border = "#333333",
                 fill = "white", cex = 1.0, subcex = 0.92) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = fill, border = border, lwd = lwd)
  if (!is.null(sub)) {
    text(x, y + h/2 - 1.7, lab, cex = cex, font = 2, adj = 0.5)
    text(x, y - 1.4, sub, cex = subcex, adj = 0.5)
  } else text(x, y, lab, cex = cex, adj = 0.5)
}
arw <- function(x1, y1, x2, y2) {
  arrows(x1, y1, x2, y2, length = 0.07, lwd = 1.1, col = "#333333")
}

# Identification
boxx(50, 121, 60, 6, "Identification", fill = "#F0F0F0", cex = 1.05)
boxx(23, 108, 42, 12, "Records identified from:",
     "PubMed 426 | Web of Science 220 | Embase 278", fill = "#EAF2F8",
     subcex = 0.82)
boxx(77, 108, 42, 12, "Records removed before screening:",
     "Duplicates removed (n = 363)", fill = "#EAF2F8")
arw(44, 108, 56, 108)
arw(23, 102, 23, 96.3); arw(77, 102, 77, 96.3)

# Screening
boxx(50, 94, 60, 6, "Screening", fill = "#F0F0F0", cex = 1.05)
boxx(23, 81, 42, 12, "Records screened", "(n = 561)", fill = "#EAF2F8")
boxx(77, 81, 42, 12, "Records excluded",
     "(n = 532), incl. 1 retraction", fill = "#EAF2F8")
arw(44, 81, 56, 81)
arw(23, 75, 23, 69.3); arw(77, 75, 77, 69.3)

# Eligibility
boxx(50, 67, 60, 6, "Eligibility", fill = "#F0F0F0", cex = 1.05)
boxx(23, 54, 42, 12, "Full texts assessed", "(n = 29)", fill = "#EAF2F8")
boxx(77, 54, 42, 14, "Full texts excluded (n = 6)",
     "no transcriptomic dataset 2\nno TBI contrast 1\nnot primary research 1\noutcome not gene-level 1",
     fill = "#EAF2F8", cex = 0.95, subcex = 0.85)
arw(44, 54, 56, 54)
arw(23, 48, 23, 42.3); arw(77, 48, 77, 42.3)

# Included
boxx(50, 40, 60, 6, "Included", fill = "#F0F0F0", cex = 1.05)
boxx(38, 30, 48, 11, "Studies included in the review", "(n = 23)",
     fill = "#D6E9F7", lwd = 1.6)
boxx(38, 16, 48, 9, "Datasets reanalysed", "(n = 6)",
     fill = "#D6E9F7", lwd = 1.6)
arw(38, 24.5, 38, 20.5)
text(50, 6, "PRISMA 2020 flow diagram. From: Page MJ, et al. BMJ 2021;372:n71.",
     cex = 0.78, col = "#555555")
dev.off()

## ---------------- Fig 2. H1 - DEG burden ----------------
pngf("Fig2", 8, 5)
op <- par(mar = c(5.5, 5.8, 2.6, 2.2)); on.exit(par(op))
d <- de[order(-de$prop_de), ]
labs <- sprintf("GSE%s\n(n=%s)", short(d$dataset_id),
                format(d$n_genes_tested, big.mark = ","))
bp <- barplot(100 * d$prop_de, names.arg = labs,
              col = ifelse(d$prop_de == 0, col_grey, COL[d$dataset_id]),
              border = NA, ylim = c(0, 64), las = 1, cex.names = 0.92,
              ylab = "", xlab = "", axes = FALSE, space = 0.6)
axis(2, at = seq(0, 60, 10), las = 1, cex.axis = 1.0)
mtext("Differentially expressed genes (% of tested genes)", side = 2,
      line = 3.9, cex = 1.05)
segments(0, 0, max(bp) + 1, 0, col = "#333333")
for (i in seq_len(nrow(d)))
  text(bp[i], 100 * d$prop_de[i] + 3.0,
       sprintf("%s\n(%.1f%%)", format(d$n_de[i], big.mark = ","), 100 * d$prop_de[i]),
       cex = 0.95, adj = 0.5)
text(max(bp) - 0.4, 60, "H1: burden is not comparable\nacross datasets (0% - 49.8%)",
     cex = 0.95, adj = 1, col = "#333333")
dev.off()

## ---------------- Fig 3. H2 / H3 - list overlap ----------------
pngf("Fig3", 9, 5)
layout(matrix(c(1, 2), 1, 2), widths = c(1.15, 1))
op <- par(mar = c(5.2, 5.2, 2.6, 2.0)); on.exit(par(op))

pp <- l2[is.finite(l2$jaccard), ]; pp <- pp[order(pp$jaccard), ]
pp$lab <- gsub("_vs_", " vs ", gsub("GSE", "", pp$pair))
bp <- barplot(pp$jaccard, horiz = TRUE, names.arg = pp$lab, col = col_main,
              border = NA, las = 1, cex.names = 0.88, xlim = c(0, 0.155),
              xlab = "", axes = FALSE)
axis(1, at = seq(0, 0.15, 0.05), labels = sprintf("%.2f", seq(0, 0.15, 0.05)), cex.axis = 0.92)
mtext("Jaccard index (full DEG lists)", side = 1, line = 2.9, cex = 1.0)
text(pp$jaccard + 0.004, bp, sprintf("%.3f", pp$jaccard), cex = 0.82, adj = 0)
mtext("A  H2: DEG-list membership overlap", side = 3, line = 0.8, cex = 1.05, adj = 0, font = 2)

ks <- sort(unique(l3$k)); kk <- l3[l3$pair %in% pp$pair, ]
plot(NA, xlim = range(ks), ylim = c(0, max(kk$overlap, na.rm = TRUE) * 1.15),
     xlab = "", ylab = "", axes = FALSE, log = "x")
axis(1, at = ks, labels = format(ks, big.mark = ","), cex.axis = 0.92)
axis(2, at = seq(0, 0.5, 0.1), labels = sprintf("%.1f", seq(0, 0.5, 0.1)), las = 1, cex.axis = 0.92)
mtext("k (top-ranked genes)", side = 1, line = 2.9, cex = 1.0)
mtext("Overlap of top-k lists", side = 2, line = 3.4, cex = 1.0)
for (pr in unique(kk$pair)) {
  s <- kk[kk$pair == pr, ]; s <- s[order(s$k), ]
  lines(s$k, s$overlap, col = "#9BB7CC", lwd = 1)
}
s2 <- kk[kk$pair == "GSE269748_vs_GSE45997", ]; s2 <- s2[order(s2$k), ]
lines(s2$k, s2$overlap, col = col_acc, lwd = 2.4)
points(s2$k, s2$overlap, col = col_acc, pch = 16, cex = 1.1)
box()
legend("topleft", legend = c("all dataset pairs", "GSE269748 vs GSE45997 (highest)"),
       col = c("#9BB7CC", col_acc), lwd = c(1, 2.4), bty = "n", cex = 0.85)
mtext("B  H3: top-k list overlap", side = 3, line = 0.8, cex = 1.05, adj = 0, font = 2)
dev.off()

## ---------------- Fig 4. H4 concordance & internal stability ----------------
pngf("Fig4", 9, 5.6)
layout(matrix(c(1, 2), 1, 2), widths = c(1.15, 1))
op <- par(mar = c(8.2, 4.6, 2.6, 1.6)); on.exit(par(op))

c4 <- l4[!is.na(l4$concordance_rate), ]; c4 <- c4[order(-c4$concordance_rate), ]
c4$lab <- gsub("_vs_", " vs ", gsub("GSE", "", c4$pair))
n2 <- c4$n_double_sig
ci <- t(mapply(wil, round(c4$concordance_rate * n2), n2))
bp <- barplot(c4$concordance_rate, names.arg = c4$lab, col = col_main, border = NA,
              ylim = c(0, 1.06), las = 2, cex.names = 0.78, axes = FALSE,
              ylab = "", xlab = "")
axis(2, at = seq(0, 1, 0.2), labels = sprintf("%.1f", seq(0, 1, 0.2)), las = 1, cex.axis = 0.92)
mtext("Concordance rate", side = 2, line = 3.2, cex = 1.0)
arrows(bp, ci[, 1], bp, ci[, 2], code = 3, angle = 90, length = 0.05,
       col = "#333333", lwd = 1.1)
text(bp, 1.02, sprintf("n=%d", n2), cex = 0.78, adj = 0.5)
abline(h = 0.5, lty = 2, col = col_grey)
mtext("A  H4: directional concordance", side = 3, line = 0.8, cex = 1.05, adj = 0, font = 2)

st <- stab
st$prop <- suppressWarnings(as.numeric(st$prop_survived)); st$prop[is.na(st$prop)] <- 0
st <- st[order(-st$prop), ]
st$lab <- paste0("GSE", short(st$dataset_id))
cols <- ifelse(st$feasible, ifelse(st$prop > 0.05, col_ok, col_acc), col_grey)
bp <- barplot(st$prop, names.arg = st$lab, col = cols, border = NA,
              ylim = c(0, 1.06), las = 2, cex.names = 0.82, axes = FALSE,
              ylab = "", xlab = "")
axis(2, at = seq(0, 1, 0.2), labels = sprintf("%.1f", seq(0, 1, 0.2)), las = 1, cex.axis = 0.92)
mtext("Proportion surviving", side = 2, line = 3.2, cex = 1.0)
text(bp, st$prop + 0.05, sprintf("%.1f%%", 100 * st$prop), cex = 0.82, adj = 0.5)
legend("topright", legend = c("feasible", "not feasible"), fill = c(col_ok, col_grey),
       border = NA, bty = "n", cex = 0.82)
mtext("B  Internal stability (SS2.7)", side = 3, line = 0.8, cex = 1.05, adj = 0, font = 2)
dev.off()

## ---------------- Fig 5. H5 by target dataset ----------------
pngf("Fig5", 8, 5)
op <- par(mar = c(5.8, 5.4, 2.8, 2.4)); on.exit(par(op))
d5 <- l5d[order(-l5d$confirm_rate), ]
labs <- sprintf("GSE%s\n(n=%d/%d)", short(d5$target_dataset),
                d5$n_confirmed, d5$n_genes_mapped)
bp <- barplot(100 * d5$confirm_rate, names.arg = labs,
              col = ifelse(d5$confirm_rate == 0, col_grey, COL[d5$target_dataset]),
              border = NA, ylim = c(0, 80), las = 1, cex.names = 0.88,
              axes = FALSE, ylab = "", xlab = "", space = 0.6)
axis(2, at = seq(0, 75, 15), las = 1, cex.axis = 1.0)
mtext("Confirmed nominations (%)", side = 2, line = 3.7, cex = 1.05)
segments(0, 0, max(bp) + 1, 0, col = "#333333")
arrows(bp, 100 * d5$ci_lo, bp, 100 * d5$ci_hi, code = 3, angle = 90,
       length = 0.05, col = "#333333", lwd = 1.2)
text(bp, 100 * pmax(d5$ci_hi, d5$confirm_rate) + 4.2,
     sprintf("%.1f%%", 100 * d5$confirm_rate), cex = 0.95, adj = 0.5)
abline(h = 25.5, lty = 2, col = col_acc, lwd = 1.2)
abline(h = c(15, 40), lty = 3, col = "#555555")
text(max(bp) - 0.3, 28.5, "pooled 25.5%", cex = 0.88, col = col_acc, adj = 1)
text(max(bp) - 0.3, 43.0, "a-priori 40%", cex = 0.8, col = "#555555", adj = 1)
text(max(bp) - 0.3, 12.0, "a-priori 15%", cex = 0.8, col = "#555555", adj = 1)
dev.off()

## ---------------- Fig 6. H5 forest plot by nominating study ----------------
pngf("Fig6", 8, 10.5)
op <- par(mar = c(4.6, 3.4, 3.2, 8.4)); on.exit(par(op))
s5 <- l5s[order(-l5s$confirm_rate, l5s$study_id), ]
n <- nrow(s5)
plot(NA, xlim = c(0, 0.62), ylim = c(0.4, n + 2.6), axes = FALSE, xlab = "", ylab = "")
axis(1, at = seq(0, 0.6, 0.1), labels = sprintf("%.0f%%", 100 * seq(0, 0.6, 0.1)), cex.axis = 0.92)
mtext("Confirmed nominations", side = 1, line = 2.7, cex = 1.0)
abline(v = seq(0, 0.6, 0.1), col = "#EEEEEE")
abline(v = 0.255, col = col_acc, lty = 2, lwd = 1.2)
for (i in seq_len(n)) {
  y <- n - i + 1.6; r <- s5[i, ]
  arrows(r$ci_lo, y, r$ci_hi, y, code = 3, angle = 90, length = 0.035,
         col = "#777777", lwd = 1.0)
  points(r$confirm_rate, y, pch = 22,
         bg = ifelse(r$confirm_rate == 0, col_grey, col_main), col = "#333333", cex = 1.5)
  mtext(r$study_id, side = 2, at = y, line = 0.8, las = 1, cex = 0.82, adj = 1)
}
segments(0, 0.9, 0.62, 0.9)
text(0.63, n + 2.1, "confirmed / mapped", cex = 0.82, adj = 0, font = 2, xpd = NA)
for (i in seq_len(n)) {
  y <- n - i + 1.6; r <- s5[i, ]
  text(0.63, y, sprintf("%d/%d", r$n_confirmed, r$n_genes_mapped), cex = 0.8, adj = 0, xpd = NA)
}
text(0.255, 0.15, "pooled 25.5%", cex = 0.85, col = col_acc, adj = 0.5, xpd = NA)
mtext("Nominating study", side = 2, line = 2.4, cex = 1.0, las = 3)
mtext("H5: confirmation by nominating study (95% Wilson CI)", side = 3,
      line = 1.0, cex = 1.05, font = 2, adj = 0)
dev.off()

## ---------------- Fig 7. Reporting practice ----------------
pngf("Fig7", 8, 5.2)
op <- par(mar = c(6.8, 4.8, 2.6, 1.8)); on.exit(par(op))
items <- c("L5: names\ngenes", "L7: reports\nthreshold", "Contrast\nstated",
           "Direction\nstated", "Subject ID\n(rows)", "Stability-elig.\n(rows)")
pct  <- c(95.7, 8.7, 95.7, 78.3, 23.3, 38.3)
frac <- c("22/23", "2/23", "22/23", "18/23", "14/60", "23/60")
bp <- barplot(pct, names.arg = items, col = ifelse(pct >= 50, col_main, col_acc),
              border = NA, ylim = c(0, 112), axes = FALSE, ylab = "", xlab = "",
              cex.names = 0.9, space = 0.7)
axis(2, at = seq(0, 100, 20), las = 1, cex.axis = 0.95)
mtext("% of units", side = 2, line = 3.2, cex = 1.0)
segments(0, 0, max(bp) + 1, 0, col = "#333333")
text(bp, pct + 9, sprintf("%.1f%%", pct), cex = 0.95, adj = 0.5)
text(bp, pct + 3.2, frac, cex = 0.78, col = "#444444", adj = 0.5)
mtext("L5 = conclusion-level gene nomination; L7 = explicit decision threshold.",
      side = 1, line = 4.6, cex = 0.78, col = "#555555")
dev.off()

cat("PDF figures written to", OUTPDF, "\n")
