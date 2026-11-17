##########################################################################
#### 02  Manuscript figures (no titles; legends in figures/Figure_legends.md)
#### Needs results/dom96/ from 01_ld_analysis.R. Each figure is written as PDF and 600-dpi PNG.
####   Main figure:          excess ALD over the permutation null, 96 domesticated lines with copy number (dom96)
####                         figures/Main/Main_ExcessALD_heatmap_dom96_{forward,reverse,paired}
####   Supplementary figure: read depth vs copy number, a6 coordinates (dom96)
####                         figures/Supplementary_Manuscript/Supp_ReadDepth_vs_CN_a6_dom96
#### Genotypes are used as called (Beagle-phased; each homolog = one haplotype); heterozygous calls are NOT set to missing.
##########################################################################

source("R/00_settings_functions.R")
dir.create("figures/Main", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/Supplementary_Manuscript", recursive = TRUE, showWarnings = FALSE)

### ---- Main figure: excess ALD heatmap ----------------------------------------------------------------------
### Excess = ALD - permutation-null mean for the direction shown (forward: ALD(j | i); reverse: ALD(i | j)); dot = FDR < 0.05
excess_figure <- function(set, dir, prefix) {
  D <- load_set(set); L <- select_loci(D)
  P <- read.csv(file.path("results", set, "Pairs_permutation.csv")); stopifnot(identical(P$Locus_i, L$Lab[P$i]))
  f <- function(x) file.path("figures", dir, paste0(prefix, x))
  eb <- max(abs(quantile(c(P$Excess_i_given_j, P$Excess_j_given_i), c(0.01, 0.99)))); brk <- seq(-eb, eb, length.out = 22)
  ex <- function(rev) { d <- if (rev) "i_given_j" else "j_given_i"
    draw_triangle(P[[paste0("Excess_", d)]], P, L, rev, excess_palette, brk, ald_key(rev, "Excess "), mark = P[[paste0("FDR_", d)]] < 0.05) }
  for (rev in c(FALSE, TRUE)) save_fig(f(paste0("ExcessALD_heatmap_", set, if (rev) "_reverse" else "_forward")), function() ex(rev))
  save_fig(f(paste0("ExcessALD_heatmap_", set, "_paired")), w = 7, h = 9.5, function() {
    par(mfrow = c(2, 1), cex = 1)                                  # mfrow would otherwise shrink all text to 0.83
    ex(FALSE); mtext("a", side = 3, adj = 0, line = -1.4, font = 2, cex = 1.4); ex(TRUE); mtext("b", side = 3, adj = 0, line = -1.4, font = 2, cex = 1.4) })
}

### ---- Supplementary figure: read depth vs copy number in a6 coordinates (2022 reads, a2 alignment, positions transposed to a6)
readdepth_figure <- function(set, dir, prefix) {
  D <- load_set(set); B <- D$Blocks
  R <- read.csv("data/readdepth_a6.csv.gz", check.names = FALSE)
  cn <- setNames(D$Lines$Copy_Number, D$Lines$Line_ID); ids <- intersect(names(cn)[!is.na(cn)], names(R))
  Dp <- as.matrix(R[, ids]); pos <- R$Pos_a6; cn <- cn[ids]
  base <- abs(pos - mean(Slocus)) > 1e6                            # depth normalised per line on SNPs > 1 Mb from the S-locus
  Rel <- sweep(Dp, 2, colSums(Dp[base, ]) / sum(base), "/")
  Br <- seq(38850000, 39350000, by = 2000); bin <- cut(pos, Br, labels = FALSE)
  Bins <- do.call(rbind, lapply(sort(unique(bin[!is.na(bin)])), function(b) { rr <- colMeans(Rel[which(bin == b), , drop = FALSE])
    data.frame(Start = Br[b], End = Br[b + 1], nSNPs = sum(bin == b, na.rm = TRUE),
               Rho_CN = suppressWarnings(cor(rr, cn, method = "spearman")), Slope_per_copy = unname(coef(lm(rr ~ cn))[2])) }))
  write.csv(Bins, file.path("results", set, "ReadDepth_vs_CN_2kb_bins.csv"), row.names = FALSE)
  ga <- t(sapply(GA2ox8, gene_pos)); zl <- c(Slocus["Start"] - 35000, Slocus["End"] + 35000)
  Bz <- B[B$BP2 >= zl[1] & B$BP1 <= zl[2] & !B$Block %in% artefact_blocks(D), ]   # artefact block (Hap119 in dom96) omitted
  draw <- function() {
    par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 1), mgp = c(2.6, 0.7, 0), cex = 1)   # ~10 pt text when 7 in wide is printed at A4 text width
    plot((Bins$Start + 1000) / 1e6, Bins$Rho_CN, type = "h", lwd = 1.5, col = ifelse(Bins$Rho_CN >= 0.5, "red3", "grey50"), ylim = c(-0.4, 1),
         xlab = "Gm13 position (Mb, Wm82.a6.v1)", ylab = expression("Spearman " * rho * " (depth, CN)"))
    rect(Slocus["Start"] / 1e6, -2, Slocus["End"] / 1e6, 2, col = adjustcolor("gold", 0.3), border = NA); abline(h = c(0, 0.5), lty = 2)
    mtext("a", side = 3, adj = -0.09, line = -0.1, font = 2, cex = 1.3)
    zb <- Bins$Start >= zl[1] & Bins$End <= zl[2]
    plot((Bins$Start[zb] + 1000) / 1e6, Bins$Slope_per_copy[zb], type = "h", lwd = 3, col = "steelblue", xlim = zl / 1e6, ylim = c(-0.45, 1.45),
         xlab = "Gm13 position (Mb, Wm82.a6.v1)", ylab = "Relative depth per copy")
    rect(Slocus["Start"] / 1e6, -2, Slocus["End"] / 1e6, 2, col = adjustcolor("gold", 0.3), border = NA)
    rect(ga[, 1] / 1e6, -2, ga[, 2] / 1e6, 2, col = adjustcolor("purple", 0.3), border = NA)
    text(rowMeans(ga) / 1e6, 1.37, rownames(ga), cex = 0.9, font = 3, col = "purple")
    abline(v = CNV_unit / 1e6, lty = 3, col = "red3")
    yb <- rep(c(-0.2, -0.32), length.out = nrow(Bz)); spanning <- Bz$BP1 < Slocus["Start"] & Bz$BP2 > Slocus["End"]
    inside <- Bz$BP1 >= Slocus["Start"] & Bz$BP1 <= Slocus["End"]
    if (nrow(Bz)) { segments(pmax(Bz$BP1, zl[1]) / 1e6, yb, pmin(Bz$BP2, zl[2]) / 1e6, yb, lwd = 4, col = ifelse(inside | spanning, "red3", "grey40"))
    text((pmax(Bz$BP1, zl[1]) + pmin(Bz$BP2, zl[2])) / 2e6, yb - 0.08, Bz$Block, cex = 0.8) }
    abline(h = 0, lty = 2); mtext("b", side = 3, adj = -0.09, line = -0.1, font = 2, cex = 1.3)
  }
  save_fig(file.path("figures", dir, paste0(prefix, "ReadDepth_vs_CN_a6_", set)), draw, w = 7, h = 6)
  cat(sprintf("%s read depth: %d lines with CN; bins with rho >= 0.5 near the S-locus: %s\n", set, length(ids),
              paste(range(Bins$Start[Bins$Rho_CN >= 0.5 & Bins$nSNPs >= 3 & Bins$Start > Slocus["Start"] - 5e4 & Bins$End < Slocus["End"] + 5e4]), collapse = "-")))
}

excess_figure("dom96", "Main", "Main_")
readdepth_figure("dom96", "Supplementary_Manuscript", "Supp_")

### Supplementary tables: the CNV row of each line set (not a figure; kept here so run_all.R writes them)
dir.create("tables", showWarnings = FALSE)
file.copy("results/dom96/CNV_vs_blocks.csv", "tables/Table_S1_CNV_vs_blocks_dom96.csv", overwrite = TRUE)
file.copy("results/all106/CNV_vs_blocks.csv", "tables/Table_S2_CNV_vs_blocks_all106.csv", overwrite = TRUE)
file.copy("results/dom99/CNV_vs_blocks.csv", "tables/Table_S3_CNV_vs_blocks_dom99.csv", overwrite = TRUE)
