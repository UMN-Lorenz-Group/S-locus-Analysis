##########################################################################
#### 05  Chance level of ALD and the chance-corrected (adjusted) ALD
#### Needs results/<set>/ from 01_ld_analysis.R (loci.csv, Pairs_permutation.csv, CNV_vs_blocks.csv).
####
#### ALD(A | B)^2 equals Goodman & Kruskal's tau(A | B). Under independence E[tau] ~ (k_B - 1) / (N - 1)
#### (k_B = alleles at the conditioning locus B, N = haplotypes; Light & Margolin 1971), so the chance level of
#### ALD is ~ sqrt((k_B - 1) / (N - 1)) and rises with the number of haplotypes at B.
#### Adjusted ALD = (ALD - null) / (1 - null), the AMI-style scaling with max(ALD) = 1; null = permutation-null mean.
####
#### Outputs (figures/Supplementary/, PDF + 600-dpi PNG):
####   <prefix>ChanceCorrection_<set>      a  permutation-null mean vs sqrt((k_B - 1) / (N - 1)), every pair, both directions
####                                       b  ALD(CNV | block) and its null vs no. of block haplotypes
####                                       c  excess and adjusted ALD(CNV | block) vs no. of block haplotypes
####   <prefix>AdjustedALD_heatmap_<set>_paired   adjusted ALD, forward (a) and reverse (b) order; dot = FDR < 0.05
####   results/<set>/Null_vs_tau_expectation.csv  the panel-a values
##########################################################################

source("R/00_settings_functions.R")

chance_figures <- function(set, dir, prefix) {
  r <- file.path("results", set); D <- load_set(set); L <- select_loci(D)
  Loc <- read.csv(file.path(r, "loci.csv")); P <- read.csv(file.path(r, "Pairs_permutation.csv"))
  P$Adj_i_given_j <- P$Excess_i_given_j / (1 - P$Null_i_given_j); P$Adj_j_given_i <- P$Excess_j_given_i / (1 - P$Null_j_given_i)
  C <- read.csv(file.path(r, "CNV_vs_blocks.csv"))
  N <- 2 * sum(!is.na(D$Lines$Copy_Number)); k <- Loc$nHaplotypes; kc <- which(L$Loci == 0)
  f <- function(x) file.path("figures", dir, paste0(prefix, x))

  ### a: every pair, both directions; the conditioning locus sets the expectation
  E <- rbind(data.frame(Pair = paste(P$Locus_i, P$Locus_j), Direction = "i_given_j", k_cond = k[P$j], Null = P$Null_i_given_j,
                        CNV = P$j == kc | P$i == kc, CNV_conditioning = P$j == kc),
             data.frame(Pair = paste(P$Locus_i, P$Locus_j), Direction = "j_given_i", k_cond = k[P$i], Null = P$Null_j_given_i,
                        CNV = P$j == kc | P$i == kc, CNV_conditioning = P$i == kc))
  E$Expected <- sqrt((E$k_cond - 1) / (N - 1)); E$Ratio <- E$Null / E$Expected
  write.csv(E, file.path(r, "Null_vs_tau_expectation.csv"), row.names = FALSE)
  rr <- cor(E$Null, E$Expected); mr <- median(E$Ratio)
  Ne <- 1 + median((E$k_cond - 1) / E$Null^2)                       # effective N implied by the null: E = sqrt((k-1)/(Ne-1))

  sp <- function(a, b) cor(a, b, method = "spearman")
  obs <- C$Null_CNV_given_block + C$Excess_CNV_given_block          # haplotype-count ALD, the one the null and excess refer to
  draw <- function() {
    par(mfrow = c(1, 3), mar = c(4, 4, 1.5, 0.8), mgp = c(2.3, 0.6, 0), cex = 0.75, las = 1)
    o <- order(E$CNV)
    plot(E$Expected[o], E$Null[o], pch = 16, cex = ifelse(E$CNV[o], 0.55, 0.35),
         col = ifelse(E$CNV[o], ifelse(E$CNV_conditioning[o], "darkorange2", "red3"), adjustcolor("grey40", 0.35)),
         xlim = c(0, max(E$Expected, E$Null)), ylim = c(0, max(E$Expected, E$Null)),
         xlab = expression(sqrt((italic(k)[B] - 1) / (italic(N) - 1))), ylab = "Permutation-null mean ALD(A | B)")
    abline(0, 1, lty = 2); abline(0, mr, col = "steelblue", lwd = 1.2)
    legend("topleft", bty = "n", cex = 0.8, pch = c(16, 16, 16, NA, NA), lty = c(NA, NA, NA, 2, 1),
           col = c("grey40", "red3", "darkorange2", "black", "steelblue"),
           legend = c("block pairs", "CNV | block", "block | CNV", "1:1", sprintf("median ratio %.2f", mr)))
    mtext(sprintf("r = %.3f, %d values", rr, nrow(E)), side = 1, line = -1.3, adj = 0.97, cex = 0.6)
    mtext("a", side = 3, adj = -0.2, line = 0.2, font = 2)

    kk <- seq(min(C$nHaplotypes), max(C$nHaplotypes), length.out = 200)
    plot(C$nHaplotypes, obs, log = "x", pch = 16, cex = 0.6, ylim = c(0, 1), xlab = "Block haplotypes",
         ylab = "ALD(CNV | block)")
    points(C$nHaplotypes, C$Null_CNV_given_block, pch = 1, cex = 0.6, col = "grey35")
    lines(kk, sqrt((kk - 1) / (N - 1)), lty = 2)
    legend("topleft", bty = "n", cex = 0.8, pch = c(16, 1, NA), lty = c(NA, NA, 2), col = c("black", "grey35", "black"),
           legend = c(sprintf("observed (rho = %.2f)", sp(C$nHaplotypes, obs)), "permutation-null mean", "tau expectation"))
    mtext("b", side = 3, adj = -0.2, line = 0.2, font = 2)

    plot(C$nHaplotypes, C$Excess_CNV_given_block, log = "x", pch = 16, cex = 0.6, col = "steelblue",
         ylim = c(min(0, C$Excess_CNV_given_block, C$Adj_CNV_given_block), 1.3 * max(C$Excess_CNV_given_block, C$Adj_CNV_given_block)), xlab = "Block haplotypes", ylab = "Chance-corrected ALD(CNV | block)")
    points(C$nHaplotypes, C$Adj_CNV_given_block, pch = 17, cex = 0.6, col = "red3"); abline(h = 0, lty = 2)
    legend("topleft", bty = "n", cex = 0.8, pch = c(16, 17), col = c("steelblue", "red3"),
           legend = c(sprintf("ALD - null (rho = %.2f)", sp(C$nHaplotypes, C$Excess_CNV_given_block)),
                      sprintf("(ALD - null) / (1 - null) (rho = %.2f)", sp(C$nHaplotypes, C$Adj_CNV_given_block))))
    mtext("c", side = 3, adj = -0.2, line = 0.2, font = 2)
  }
  save_fig(f(paste0("ChanceCorrection_", set)), draw, w = 7.5, h = 2.9)

  ### Adjusted-ALD heatmaps, same layout as the excess heatmaps
  eb <- max(abs(quantile(c(P$Adj_i_given_j, P$Adj_j_given_i), c(0.01, 0.99)))); brk <- seq(-eb, eb, length.out = 22)
  ad <- function(rev) { d <- if (rev) "i_given_j" else "j_given_i"
    draw_triangle(P[[paste0("Adj_", d)]], P, L, rev, excess_palette, brk, ald_key(rev, "Adjusted "), mark = P[[paste0("FDR_", d)]] < 0.05) }
  save_fig(f(paste0("AdjustedALD_heatmap_", set, "_paired")), w = 7, h = 9.5, function() {
    par(mfrow = c(2, 1)); ad(FALSE); mtext("a", side = 3, adj = 0, line = -1.2, font = 2); ad(TRUE); mtext("b", side = 3, adj = 0, line = -1.2, font = 2) })

  cat(sprintf("\n== %s chance level: N = %d haplotypes; null vs sqrt((k_B-1)/(N-1)) r = %.3f, median ratio %.3f (implied effective N %.0f)\n",
              set, N, rr, mr, Ne))
  cat(sprintf("   CNV | block, Spearman with no. of haplotypes: observed %.2f, excess %.2f, adjusted %.2f; excess vs adjusted rho %.3f\n",
              sp(C$nHaplotypes, obs), sp(C$nHaplotypes, C$Excess_CNV_given_block), sp(C$nHaplotypes, C$Adj_CNV_given_block),
              sp(C$Excess_CNV_given_block, C$Adj_CNV_given_block)))
  cat(sprintf("   block | CNV: excess vs adjusted rho %.3f; adjusted range %.3f to %.3f\n",
              sp(C$Excess_block_given_CNV, C$Adj_block_given_CNV), min(C$Adj_block_given_CNV), max(C$Adj_block_given_CNV)))
  invisible(E)
}

chance_figures("dom96", "Supplementary", "Supp_")
chance_figures("all106", "Supplementary", "Supp_Soja_")
chance_figures("dom99", "Supplementary", "Supp_Blocks99_")
