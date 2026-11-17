##########################################################################
#### 03  Summary statistics for the Results text: S-locus CNV vs haploblocks
#### Needs results/<set>/CNV_vs_blocks.csv from 01_ld_analysis.R.
#### Output: tables/Summary_CNV_vs_blocks_<set>.csv (median, mean, SD, IQR, range per measure)
##########################################################################

source("R/00_settings_functions.R")
dir.create("tables", showWarnings = FALSE)

summarise_cnv <- function(set) {
  C <- read.csv(file.path("results", set, "CNV_vs_blocks.csv"))
  Measures <- c(ALD_CNV_given_block    = "ALD(CNV | block), haplo.em",
                ALD_block_given_CNV    = "ALD(block | CNV), haplo.em",
                Null_CNV_given_block   = "Permutation-null mean, ALD(CNV | block)",
                Null_block_given_CNV   = "Permutation-null mean, ALD(block | CNV)",
                Excess_CNV_given_block = "Excess ALD(CNV | block)",
                Excess_block_given_CNV = "Excess ALD(block | CNV)",
                Adj_CNV_given_block    = "Adjusted ALD(CNV | block) = (ALD - null) / (1 - null)",
                Adj_block_given_CNV    = "Adjusted ALD(block | CNV) = (ALD - null) / (1 - null)",
                AMI                    = "Adjusted mutual information")
  S <- do.call(rbind, lapply(names(Measures), function(m) { v <- C[[m]]
    data.frame(Measure = Measures[[m]], n_blocks = length(v), Median = median(v), Mean = mean(v), SD = sd(v),
               Q1 = unname(quantile(v, 0.25)), Q3 = unname(quantile(v, 0.75)), Min = min(v), Max = max(v),
               Block_at_min = C$Locus[which.min(v)], Block_at_max = C$Locus[which.max(v)]) }))
  write.csv(S, file.path("tables", paste0("Summary_CNV_vs_blocks_", set, ".csv")), row.names = FALSE)

  cat(sprintf("\n== %s: %d blocks vs the S-locus CNV (%d upstream, %d downstream)\n", Sets[[set]]$Label, nrow(C),
              sum(C$Side == "upstream"), sum(C$Side == "downstream")))
  print(S[, c("Measure", "Median", "Mean", "SD", "Min", "Max", "Block_at_max")], row.names = FALSE, digits = 3)
  sp <- function(a, b) cor(a, b, method = "spearman")
  cat(sprintf("ALD(CNV | block) > ALD(block | CNV) in %d of %d blocks\n", sum(C$ALD_CNV_given_block > C$ALD_block_given_CNV), nrow(C)))
  cat(sprintf("Blocks with ALD(CNV | block) >= 0.5: %d (median %g haplotypes, vs %g for the rest)\n", sum(C$ALD_CNV_given_block >= 0.5),
              median(C$nHaplotypes[C$ALD_CNV_given_block >= 0.5]), median(C$nHaplotypes[C$ALD_CNV_given_block < 0.5])))
  cat(sprintf("Spearman with no. of block haplotypes: ALD(CNV | block) %.2f -> excess %.2f; ALD(block | CNV) %.2f -> excess %.2f\n",
              sp(C$nHaplotypes, C$ALD_CNV_given_block), sp(C$nHaplotypes, C$Excess_CNV_given_block),
              sp(C$nHaplotypes, C$ALD_block_given_CNV), sp(C$nHaplotypes, C$Excess_block_given_CNV)))
  cat(sprintf("  adjusted (ALD - null) / (1 - null): ALD(CNV | block) %.2f; ALD(block | CNV) %.2f\n",
              sp(C$nHaplotypes, C$Adj_CNV_given_block), sp(C$nHaplotypes, C$Adj_block_given_CNV)))
  cat(sprintf("FDR < 0.05: ALD(CNV | block) %d of %d; ALD(block | CNV) %d of %d\n",
              sum(C$FDR_CNV_given_block < 0.05), nrow(C), sum(C$FDR_block_given_CNV < 0.05), nrow(C)))
  cat(sprintf("Haplotypes per block: median %g, range %g-%g\n", median(C$nHaplotypes), min(C$nHaplotypes), max(C$nHaplotypes)))
  invisible(S)
}

for (s in names(Sets)) summarise_cnv(s)
