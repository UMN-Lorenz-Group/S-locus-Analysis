##########################################################################
#### 01  Pairwise asymmetric LD, PH13 region to the S-locus CNV, for each line set
#### Results (results/<set>/), reused when present (delete a file to recompute it):
####   loci.csv                     loci in genomic order (PLINK block, span, SNPs, haplotypes)
####   ALD_matrix_haplo.em.csv      [i, j] = ALD(locus j | locus i), [j, i] = ALD(locus i | locus j); i upstream of j
####                                (dom99: ~10 min; all106: longer)
####   Pairs_permutation.csv        every pair: ALD both directions (haplotype counts), permutation-null mean,
####                                excess, P, BH-FDR, AMI (dom99: ~5 min; all106: ~20 min)
####   CNV_vs_blocks.csv            the CNV row: ALD(CNV | block) and ALD(block | CNV) with null, excess, adjusted
####                                ((ALD - null) / (1 - null)), P, FDR; AMI
##########################################################################

source("R/00_settings_functions.R")

run_set <- function(set) {
  out <- file.path("results", set); 
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  D <- load_set(set); 
  L <- select_loci(D)
  keep <- which(!is.na(D$Lines$Copy_Number))                     # lines with a copy-number estimate
  
  G <- locus_genotypes(D, L, keep); nI <- length(keep); n <- length(L$Loci); kc <- which(L$Loci == 0)
  B <- D$Blocks
  Loc <- data.frame(Order = seq_len(n), Locus = L$Lab, Block = ifelse(L$Loci == 0, NA, L$Loci),
                    Start = ifelse(L$Loci == 0, Slocus["Start"], B$BP1[pmax(L$Loci, 1)]),
                    End = ifelse(L$Loci == 0, Slocus["End"], B$BP2[pmax(L$Loci, 1)]),
                    nSNPs = ifelse(L$Loci == 0, NA, B$NSNPS[pmax(L$Loci, 1)]), nHaplotypes = vapply(G, max, 0))
 
  cat(sprintf("\n== %s: %d lines (%d with CN), %d blocks; heatmap %s to %s, %d loci; PH13 block %d; S-locus: %s\n",
              Sets[[set]]$Label, nrow(D$Geno), nI, nrow(B), L$Lab[1], L$Lab[n], n, L$PH13_Block,
              if (length(L$Spanning)) sprintf("Hap%s spans it (kept next to the CNV)", paste(L$Spanning, collapse = ","))
              else if (length(L$Inside)) sprintf("Hap%s start inside (replaced by the CNV)", paste(range(L$Inside), collapse = "-"))
              else sprintf("no block overlaps it (CNV placed between Hap%d and Hap%d)", L$Loci[which(L$Loci == 0) - 1], L$Loci[which(L$Loci == 0) + 1])))

  fM <- file.path(out, "ALD_matrix_haplo.em.csv")
  if (file.exists(fM)){ M <- as.matrix(read.csv(fM, row.names = 1, check.names = FALSE)); stopifnot(identical(rownames(M), L$Lab))}else{
         t0 <- Sys.time();
		 M <- pairwise_ald(G, L$Lab); 
		 write.csv(M, fM)
         cat(sprintf("haplo.em pairs: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }

  fP <- file.path(out, "Pairs_permutation.csv")
  if (file.exists(fP)) { P <- read.csv(fP); stopifnot(identical(P$Locus_i, L$Lab[P$i])) } else{
         t0 <- Sys.time(); 
         P <- pair_permutation(G, L$Lab, nI); 
         write.csv(P, fP, row.names = FALSE)
         cat(sprintf("permutations: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))) 
  }
  
  ### Chance-corrected ALD, scaled like AMI with max(ALD) = 1: (ALD - null) / (1 - null) = excess / (1 - null)
  P$Adj_i_given_j <- P$Excess_i_given_j / (1 - P$Null_i_given_j);
  P$Adj_j_given_i <- P$Excess_j_given_i / (1 - P$Null_j_given_i)
  
  d <- c(abs(M[cbind(P$i, P$j)] - P$ALD_j_given_i), abs(M[cbind(P$j, P$i)] - P$ALD_i_given_j))
  cat(sprintf("haplo.em vs haplotype-count ALD: |diff| median %.4f, 95th pct %.4f\n", median(d), quantile(d, 0.95)))
  
  write.csv(Loc, file.path(out, "loci.csv"), row.names = FALSE)

  ### CNV row: [1] ALD(CNV | block), [2] ALD(block | CNV)
  C <- do.call(rbind, lapply(setdiff(seq_len(n), kc), function(b) {
    
	r <- P[(P$i == b & P$j == kc) | (P$i == kc & P$j == b), ]; 
	up <- b < kc
    g <- function(s) if (up) c(r[[paste0(s, "_j_given_i")]], r[[paste0(s, "_i_given_j")]]) else c(r[[paste0(s, "_i_given_j")]], r[[paste0(s, "_j_given_i")]])
    A <- g("ALD"); 
	N <- g("Null"); 
	X <- g("Excess"); 
	Aj <- g("Adj"); 
	Pv <- g("P"); 
	Fd <- g("FDR")
    data.frame(Locus = L$Lab[b], Start = Loc$Start[b], End = Loc$End[b], Side = if (up) "upstream" else "downstream", nHaplotypes = Loc$nHaplotypes[b],
               ALD_CNV_given_block = M[b, kc], ALD_block_given_CNV = M[kc, b],
               Null_CNV_given_block = N[1], Excess_CNV_given_block = X[1], Adj_CNV_given_block = Aj[1], P_CNV_given_block = Pv[1], FDR_CNV_given_block = Fd[1],
               Null_block_given_CNV = N[2], Excess_block_given_CNV = X[2], Adj_block_given_CNV = Aj[2], P_block_given_CNV = Pv[2], FDR_block_given_CNV = Fd[2], AMI = r$AMI)
  }))
  
  write.csv(C, file.path(out, "CNV_vs_blocks.csv"), row.names = FALSE)
  
  
  sp <- function(a, b) cor(a, b, method = "spearman")
  cat(sprintf("CNV vs %d blocks, Spearman with no. of haplotypes: ALD(CNV|block) %.2f, excess %.2f; ALD(block|CNV) %.2f, excess %.2f\n",
              nrow(C), sp(C$nHaplotypes, C$ALD_CNV_given_block), sp(C$nHaplotypes, C$Excess_CNV_given_block),
              sp(C$nHaplotypes, C$ALD_block_given_CNV), sp(C$nHaplotypes, C$Excess_block_given_CNV)))
  cat(sprintf("  adjusted (ALD - null) / (1 - null): CNV|block %.2f, block|CNV %.2f\n",
              sp(C$nHaplotypes, C$Adj_CNV_given_block), sp(C$nHaplotypes, C$Adj_block_given_CNV)))
  cat(sprintf("FDR < 0.05: CNV|block %d, block|CNV %d of %d | all pairs: %d pairs, FDR < 0.05 forward %d, reverse %d\n",
              sum(C$FDR_CNV_given_block < 0.05), sum(C$FDR_block_given_CNV < 0.05), nrow(C), nrow(P),
              sum(P$FDR_j_given_i < 0.05), sum(P$FDR_i_given_j < 0.05)))
  invisible(list(D = D, L = L, M = M, P = P, C = C))
}

Res <- lapply(setNames(names(Sets), names(Sets)), run_set) 

Res <- run_set(Sets[[1]])
SetNames <- setNames(names(Sets), names(Sets))

