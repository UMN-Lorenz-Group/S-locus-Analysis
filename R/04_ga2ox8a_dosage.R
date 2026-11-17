##########################################################################
#### 04  GA2ox8A gene dosage: read depth and allele ratios track S-locus copy number
#### (pipeline step 14, panels C and D, moved here in a6 coordinates)
#### a) mean relative read depth over the SNPs in GA2ox8A vs copy number, per line
#### b) ALT-read fraction at heterozygous calls in GA2ox8A vs copy number; a true heterozygote gives ~0.5,
####    a variant carried by one of CN paralogous copies gives ~1/CN
#### Data: data/readdepth_a6.csv.gz, data/ga2ox8a_calls_a6.csv (2022 calls, a2 alignment, positions transposed to a6)
#### Outputs: figures/Supplementary/Supp_[Soja_]GA2ox8A_Dosage_<set>.pdf/.png,
####          results/<set>/GA2ox8A_depth_by_line.csv, GA2ox8A_HetCall_ALTfraction_by_CN.csv
##########################################################################

source("R/00_settings_functions.R")

ga2ox8a_dosage_figure <- function(set, dir, prefix, min_depth = 4) {
  D <- load_set(set); g <- gene_pos(GA2ox8["GA2ox8A"])
  cn <- setNames(D$Lines$Copy_Number, D$Lines$Line_ID); cn <- cn[!is.na(cn)]

  ### a) relative depth (per line: depth / mean depth at SNPs > 1 Mb from the S-locus), averaged over the gene's SNPs
  R <- read.csv("data/readdepth_a6.csv.gz", check.names = FALSE)
  ids <- intersect(names(cn), names(R)); Dp <- as.matrix(R[, ids]); pos <- R$Pos_a6; cn <- cn[ids]
  base <- abs(pos - mean(Slocus)) > 1e6
  Rel <- sweep(Dp, 2, colSums(Dp[base, ]) / sum(base), "/")
  ig <- pos >= g[1] & pos <= g[2]
  dep <- colMeans(Rel[ig, , drop = FALSE])
  fit <- lm(dep ~ cn); rho <- cor(dep, cn, method = "spearman")
  write.csv(data.frame(Line_ID = ids, Copy_Number = cn, RelDepth_GA2ox8A = round(dep, 3)),
            file.path("results", set, "GA2ox8A_depth_by_line.csv"), row.names = FALSE)

  ### b) ALT-read fraction at het calls (depth >= min_depth), by rounded copy number (6+ pooled)
  C <- read.csv("data/ga2ox8a_calls_a6.csv")
  C <- C[C$Line_ID %in% ids & C$GT %in% c("0/1", "1/0") & C$Depth >= min_depth, ]
  C$Frac <- C$ALT_reads / C$Depth; C$CNc <- pmin(round(cn[C$Line_ID]), 6)
  cls <- sort(unique(C$CNc)); xl <- ifelse(cls == 6, "6+", cls)
  V <- data.frame(CN = xl, nHetCalls = as.vector(table(factor(C$CNc, cls))),
                  nLines = vapply(cls, function(x) length(unique(C$Line_ID[C$CNc == x])), 0L),
                  Median_ALT_frac = vapply(cls, function(x) median(C$Frac[C$CNc == x]), 0),
                  Expected_one_copy_differs = ifelse(cls == 1, NA, round(1 / cls, 3)))
  write.csv(V, file.path("results", set, "GA2ox8A_HetCall_ALTfraction_by_CN.csv"), row.names = FALSE)

  draw <- function() {
    par(mfrow = c(1, 2), mar = c(4, 4.2, 1.2, 0.8), mgp = c(2.3, 0.7, 0), cex = 0.8)
    plot(cn, dep, pch = 19, col = adjustcolor("purple", 0.6), xlab = "S-locus copy number (estimate)",
         ylab = expression("Relative read depth in " * italic(GA2ox8A)), xlim = c(0, max(cn) + 0.5), ylim = c(0, max(dep) * 1.05))
    abline(fit, lty = 2); abline(0, 1, col = "grey60", lty = 3)
    legend("topleft", c(sprintf("fit: slope %.2f per copy", coef(fit)[2]), "1:1 (depth = copy number)"), lty = c(2, 3),
           col = c("black", "grey60"), bty = "n", cex = 0.85, inset = c(0, 0.08))
    mtext(bquote(rho == .(sprintf("%.2f", rho)) * ",  n = " * .(length(ids)) * " lines"), side = 3, line = -1.2, adj = 0.04, cex = 0.8)
    mtext("a", side = 3, adj = -0.2, line = 0, font = 2)

    xs <- match(C$CNc, cls) + runif(nrow(C), -0.22, 0.22)
    plot(xs, C$Frac, pch = 16, cex = 0.5, col = adjustcolor("grey40", 0.35), xlim = c(0.5, length(cls) + 0.5), ylim = c(0, 1),
         xaxt = "n", xlab = "S-locus copy number (rounded)", ylab = "ALT-read fraction at het calls")
    axis(1, seq_along(cls), xl)
    abline(h = 0.5, lty = 2)
    segments(seq_along(cls) - 0.3, V$Median_ALT_frac, seq_along(cls) + 0.3, V$Median_ALT_frac, lwd = 3)
    points(seq_along(cls), V$Expected_one_copy_differs, pch = 4, col = "red3", cex = 1.4, lwd = 2)
    text(seq_along(cls), 0.98, paste0("n=", V$nHetCalls), cex = 0.7)
    legend("topright", c("median", "1/CN (variant in one copy)", "0.5 (true heterozygote)"), lty = c(1, NA, 2), lwd = c(3, 2, 1),
           pch = c(NA, 4, NA), col = c("black", "red3", "black"), bty = "n", cex = 0.8, inset = c(0, 0.06))
    mtext("b", side = 3, adj = -0.2, line = 0, font = 2)
  }
  save_fig(file.path("figures", dir, paste0(prefix, "GA2ox8A_Dosage_", set)), draw, w = 7, h = 3.4)
  cat(sprintf("%s GA2ox8A: %d SNPs with depth, %d lines; depth vs CN rho %.2f, slope %.2f per copy (R2 %.2f)\n",
              set, sum(ig), length(ids), rho, coef(fit)[2], summary(fit)$r.squared))
  print(V, row.names = FALSE)
}

set.seed(Seed)                                           # jitter only
ga2ox8a_dosage_figure("dom96", "Supplementary", "Supp_")
ga2ox8a_dosage_figure("all106", "Supplementary", "Supp_Soja_")
