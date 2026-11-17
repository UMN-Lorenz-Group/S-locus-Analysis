##########################################################################
#### S-locus CNV / haploblock asymmetric LD on Gm13, Wm82.a6.v1 coordinates
#### Settings and functions. Sourced by 01_ld_analysis.R and 02_figures.R.
#### Working directory = the Manuscript folder (open Slocus_Manuscript.Rproj in RStudio).
##########################################################################

suppressPackageStartupMessages({
  library(haplo.stats)   # 1.9.9.1  haplotype frequencies by EM (block-block pairs)
  library(asymLD)        # 0.1      compute.ALD (Thomson & Single 2014)
  library(aricode)       # 1.1.0    adjusted mutual information
  library(LDheatmap)     # 1.0.6    triangular heatmaps
  library(grid)
})
stopifnot("Run from the Manuscript folder (open Slocus_Manuscript.Rproj)" = file.exists("Slocus_Manuscript.Rproj"))

### ---- Settings -------------------------------------------------------------

Slocus <- c(Start = 39049348, End = 39117340)  # a2 S-locus 38,796,000-38,864,000 lifted to a6 (XMap offsets of the end SNPs) 
CNV_unit <- c(Start = 39046000, End = 39120000) # read-depth CNV unit (2-kb bins where depth tracks copy number)
PH13_Gene <- "Glyma.13G276700"
GA2ox8 <- c(GA2ox8A = "Glyma.13G287600", GA2ox8B = "Glyma.13G288000")

nPerm <- 1000                                  # permutations per pair (lines of one locus shuffled)
Seed  <- 20260928
### Line set -> PLINK 1.9 haploblocks (--blocks no-pheno-req --blocks-max-kb 500). Lines without CN are never analysed.
###   dom96  PRIMARY: blocks estimated on the 96 analysed lines, SNPs of the read-depth CNV unit removed and the
###          upstream / downstream segments run separately (no block can cross the CNV). Identical to the whole-panel
###          96-line run (S7_dom96), because no CNV-unit SNP entered a block there.
###   dom99  sensitivity: blocks on 99 domesticated lines incl. the 3 without CN (one 237-kb block spans the S-locus)
###   all106 sensitivity: 106 lines incl. 7 G. soja
Sets <- list(
  dom96  = list(Label = "96 domesticated lines with CN", Drop_Soja = TRUE, Blocks = "data/blocks/S8_dom96_noCNV_max500kb.blocks.det"),
  dom99  = list(Label = "99 domesticated lines (blocks incl. 3 lines without CN)", Drop_Soja = TRUE, Blocks = "data/blocks/S2_dom99_max500kb.blocks.det"),
  all106 = list(Label = "106 lines incl. 7 G. soja", Drop_Soja = FALSE, Blocks = "data/blocks/S1_all106_max500kb.blocks.det"))

### Heatmap span, fixed in base pairs for every line set: 38,008,447-39,820,586 (Gm13, a6). Derived from the 2022
### heatmap rule (6th block upstream of the PH13 block to the 19th block downstream of the S-locus) applied to the
### 99-line block set; it equals the start of 96-line block Hap75 and the end of Hap137.
Span <- c(38008447, 39820586)

### ---- Inputs -----------------------------------------------------------------
load_set <- function(set) {
  D <- readRDS("data/Slocus_a6_genotypes.rds")
  if (Sets[[set]]$Drop_Soja) { k <- !D$Lines$Soja; D$Geno <- D$Geno[k, , drop = FALSE]; D$Lines <- D$Lines[k, ] }
  B <- read.table(Sets[[set]]$Blocks, header = TRUE, stringsAsFactors = FALSE)
  names(B) <- c("CHR", "BP1", "BP2", "KB", "NSNPS", "SNPS"); B$Block <- seq_len(nrow(B))
  D$Blocks <- B; D$BSNPs <- lapply(strsplit(B$SNPS, "\\|"), as.integer)   # SNP indices = columns of Geno
  D$Set <- set; D
}
Genes <- read.csv("data/genes_Gm13_a6.csv")
gene_pos <- function(id) unlist(Genes[Genes$Gene == id, c("Start", "End")])

### ---- Loci: blocks in the span, CNV in place of the blocks that start inside the S-locus ----------
### dom96: no block overlaps the CNV unit, so the CNV simply sits between its flanking blocks.
### Sensitivity sets: a block that starts inside the S-locus is replaced by the CNV (2021-22 definition); a block
### spanning the whole S-locus (dom99 Hap125) is kept next to the CNV.
heatmap_span <- function() Span

select_loci <- function(D, span = heatmap_span()) {
  B <- D$Blocks
  inside <- B$BP1 >= Slocus["Start"] & B$BP1 <= Slocus["End"]
  keep <- which(B$BP1 >= span[1] & B$BP2 <= span[2] & !inside)
  loci <- c(keep[B$BP1[keep] < Slocus["Start"]], 0L, keep[B$BP1[keep] > Slocus["End"]])   # 0 = CNV
  ph <- gene_pos(PH13_Gene); kph <- which(B$BP1 <= ph[2] & B$BP2 >= ph[1])[1]
  spanning <- which(B$BP1 < Slocus["Start"] & B$BP2 > Slocus["End"])
  lab <- ifelse(loci == 0, "S-locus CNV", paste0("Hap", loci))
  lab[loci == kph] <- paste0("Hap", kph, "-PH13")
  mid <- ifelse(loci == 0, mean(Slocus), (B$BP1[pmax(loci, 1)] + B$BP2[pmax(loci, 1)]) / 2)
  list(Loci = loci, Lab = lab, Mid = mid, PH13_Block = kph, Spanning = spanning, Inside = which(inside))
}

### ---- Artefact block (read-depth figure only) ---------------------------------------------
### Hap119 (dom96; 39,149,128-39,149,590) is a fixed duplication collapsed onto one locus: its 6 SNPs are heterozygous in
### ~99% of these inbred lines (72 of 96 share one genotype string), and read depth there is ~1.7x the genome mean in every
### line, independent of S-locus copy number. Rule: median per-SNP heterozygosity > 0.5 in the lines with CN (only this
### block, in every line set; next highest 0.44). It is omitted from the block track of the read-depth figure and kept in
### the LD analysis and heatmaps.

Max_block_het <- 0.5
artefact_blocks <- function(D) {
  g <- D$Geno[!is.na(D$Lines$Copy_Number), , drop = FALSE]; h <- colMeans(substr(g, 1, 1) != substr(g, 2, 2))
  which(vapply(D$BSNPs, function(s) median(h[s]) > Max_block_het, TRUE))
}

### ---- Copy number and haplotypes -------------------------------------------------
cn3 <- function(cn) cut(round(cn), c(-Inf, 2, 4, Inf), labels = FALSE)   # <=2, 3-4, >=5 copies (few lines have >=5; the 2022 scripts used 5 classes)

### Block haplotypes: each homolog of the two-letter genotypes read as one haplotype (inbred lines, Beagle-phased
### panel), coded as integer alleles; returns lines x 2

hap_alleles <- function(Geno, snps) {
  g <- Geno[, snps, drop = FALSE]
  s1 <- apply(matrix(substr(g, 1, 1), nrow(g)), 1, paste, collapse = "")
  s2 <- apply(matrix(substr(g, 2, 2), nrow(g)), 1, paste, collapse = "")
  matrix(as.integer(factor(c(s1, s2))), ncol = 2)
}
locus_genotypes <- function(D, L, keep) {
  cn <- cn3(D$Lines$Copy_Number[keep])
  lapply(L$Loci, function(l) if (l == 0) cbind(cn, cn) else hap_alleles(D$Geno[keep, , drop = FALSE], D$BSNPs[[l]]))
}

### ---- Asymmetric LD ----------------------------------------------------------------
### ALD(A | B) = sqrt((F_A|B - F_A) / (1 - F_A)), F_A = sum p_i^2, F_A|B = sum_ij p_ij^2 / q_j (Thomson & Single 2014)
### From two-locus haplotype counts (used for the permutation null):

ald_counts <- function(a, b) {
  na <- max(a); nb <- max(b)
  p  <- matrix(tabulate((b - 1L) * na + a, na * nb), na, nb) / length(a)
  pa <- rowSums(p); pb <- colSums(p); Fa <- sum(pa^2); Fb <- sum(pb^2)
  Fab <- sum(sweep(p^2, 2, pmax(pb, 1e-300), "/")); Fba <- sum(sweep(p^2, 1, pmax(pa, 1e-300), "/"))
  c(A_given_B = sqrt(max(0, (Fab - Fa) / (1 - Fa))), B_given_A = sqrt(max(0, (Fba - Fb) / (1 - Fb))))
}

### From haplo.em two-locus haplotype frequencies + asymLD::compute.ALD (heatmap values); fixed seed, 30 starts
ald_em <- function(G1, G2, labels, seed = Seed, n.try = 30) {
  set.seed(seed)
  em <- haplo.em(geno = data.frame(G1[, 1], G1[, 2], G2[, 1], G2[, 2]), locus.label = labels,
                 control = haplo.em.control(n.try = n.try))
  hf <- cbind(em$haplotype, em$hap.prob); names(hf) <- c("allele1", "allele2", "haplo.freq")
  hf$locus1 <- labels[1]; hf$locus2 <- labels[2]
  capture.output(a <- compute.ALD(hf)); a
}

### M[i, j] = ALD(locus j | locus i) and M[j, i] = ALD(locus i | locus j); i earlier than j in genomic order.
### The upper triangle (drawn by LDheatmap) is ALD(later | earlier), as in the 2022 heatmap (ALD_HapTable2_1).
pairwise_ald <- function(G, lab) {
  n <- length(G); M <- diag(1, n)
  for (i in 1:(n - 1)) for (j in (i + 1):n) { a <- ald_em(G[[i]], G[[j]], lab[c(i, j)]); M[i, j] <- a$ALD.2.1; M[j, i] <- a$ALD.1.2 }
  dimnames(M) <- list(lab, lab); M
}
### Both directions with a permutation null (one shared set of shuffles; both homologs of a line move together)
pair_permutation <- function(G, lab, nI) {
  LC <- lapply(G, function(H) as.integer(factor(paste(pmin(H[, 1], H[, 2]), pmax(H[, 1], H[, 2])))))   # line genotype codes
  Pairs <- which(upper.tri(matrix(0, length(G), length(G))), arr.ind = TRUE); Pairs <- Pairs[order(Pairs[, 1], Pairs[, 2]), ]
  set.seed(Seed)
  
  S <- t(apply(Pairs, 1, function(p) {
    i <- p[1]; j <- p[2]; a <- c(G[[i]][, 1], G[[i]][, 2]); Hb <- G[[j]]
    obs <- ald_counts(a, c(Hb[, 1], Hb[, 2]))                    # [1] ALD(i | j), [2] ALD(j | i)
    null <- vapply(seq_len(nPerm), function(r) { o <- sample(nI); ald_counts(a, c(Hb[o, 1], Hb[o, 2])) }, numeric(2))
    
	c(ALD_i_given_j = obs[[1]], ALD_j_given_i = obs[[2]], Null_i_given_j = mean(null[1, ]), Null_j_given_i = mean(null[2, ]),
      Excess_i_given_j = obs[[1]] - mean(null[1, ]), Excess_j_given_i = obs[[2]] - mean(null[2, ]),
      P_i_given_j = (1 + sum(null[1, ] >= obs[[1]])) / (nPerm + 1), P_j_given_i = (1 + sum(null[2, ] >= obs[[2]])) / (nPerm + 1),
      AMI = AMI(LC[[i]], LC[[j]]))
  }))
  Out <- data.frame(i = Pairs[, 1], j = Pairs[, 2], Locus_i = lab[Pairs[, 1]], Locus_j = lab[Pairs[, 2]], S)
  Out$FDR_i_given_j <- p.adjust(Out$P_i_given_j, "BH"); Out$FDR_j_given_i <- p.adjust(Out$P_j_given_i, "BH")
  Out
}

### ---- Figures ---------------------------------------------------------------------
### Axis labels: block names; Mb (locus midpoint) only on the two end loci, the PH13 block and the CNV
axis_labels <- function(L) {
  n <- length(L$Loci); is_mb <- seq_len(n) %in% c(1, n) | L$Loci == 0 | L$Loci == L$PH13_Block
  txt <- ifelse(is_mb, sprintf("%s (%.2f Mb)", L$Lab, L$Mid / 1e6), L$Lab)
  tick <- which(L$Loci %in% pretty(L$Loci[L$Loci > 0], 8))                 # a few evenly spaced block names in between,
  tick <- tick[vapply(tick, function(k) min(abs(k - which(is_mb))), 0) > max(2, n %/% 25)]   # none crowding a Mb-labelled locus
  shown <- is_mb | seq_len(n) %in% tick
  list(Txt = txt, Shown = shown)
}
ald_palette <- colorRampPalette(rev(c("yellow", "orange", "red")), space = "rgb")(20)   # 2022 heatmap colours

### Colour-key labels that name the direction shown (numbered as in asymLD ALD.1.2/ALD.2.1 and PyPop W12/W21).
### Defined in the figure legend: ALD1 (forward order) = ALD of the downstream locus given the upstream locus;
### ALD2 (reverse order) = ALD of the upstream locus given the downstream locus.
ald_key <- function(rev, prefix = "") paste0(prefix, if (rev) "ALD2" else "ALD1")

### Colour key drawn with grid in the current viewport (npc units), close to the heatmap
grid_key <- function(pal, x = c(0.56, 0.80), y = 0.30, h = 0.03, title = "", at = seq(0, 1, 0.2), rng = c(0, 1), fs = 7) {
  xs <- seq(x[1], x[2], length.out = length(pal) + 1)
  grid.rect(head(xs, -1), y, diff(xs), h, just = c("left", "bottom"), gp = gpar(fill = pal, col = NA))
  grid.rect(x[1], y, diff(x), h, just = c("left", "bottom"), gp = gpar(fill = NA, col = "black", lwd = 0.5))
  xa <- x[1] + (at - rng[1]) / diff(rng) * diff(x)
  grid.segments(xa, y, xa, y - 0.008, gp = gpar(lwd = 0.5))
  grid.text(format(at), xa, y - 0.012, just = c("centre", "top"), gp = gpar(fontsize = fs - 1))
  grid.text(title, mean(x), y + h + 0.01, just = c("centre", "bottom"), gp = gpar(fontsize = fs + 1))
}

### Raw ALD heatmap (LDheatmap, own colour key); rev = TRUE draws the loci in reverse genomic order
draw_ald_heatmap <- function(M, L, rev = FALSE, label_fs = 8, key = list()) {
  n <- nrow(M); o <- if (rev) n:1 else seq_len(n); A <- axis_labels(L)
  P <- pmin(pmax(M[o, o], 0), 1); dimnames(P) <- list(A$Txt[o], A$Txt[o])    # round-off above 1 would render white
  nm <- if (rev) "ldheatmapReverse" else "ldheatmapForward"     # own grob name per panel
  LDheatmap(P, add.map = FALSE, add.key = FALSE, color = ald_palette, SNP.name = A$Txt[o][A$Shown[o]], flip = FALSE,
            newpage = FALSE, title = "", name = nm)
  grid.edit("SNPnames", gp = gpar(col = "black", fontsize = label_fs), redraw = TRUE, grep = TRUE, global = TRUE, strict = FALSE)
  do.call(grid_key, c(list(pal = rev(ald_palette), title = ald_key(rev)), key))   # LDheatmap maps the first colour to ALD = 1
}

### Triangle heatmap of a per-pair statistic (same geometry as LDheatmap); CNV row/column outlined; dots = mark
### Font sizes (cex, relative to 12 pt) chosen so that, at 7 in wide printed at A4 text width (~6.3 in), locus labels are
### ~10 pt and key labels ~9.5 pt
draw_triangle <- function(v, Pairs, L, rev, pal, brk, key_lab, mark = NULL, lab_cex = 1, key_cex = 0.95) {
  n <- length(L$Loci); pos <- if (rev) n + 1 - seq_len(n) else seq_len(n); A <- axis_labels(L); kc <- which(L$Loci == 0)
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot(NA, xlim = c(0, n * 1.07), ylim = c(-n * 0.2 * (1 + lab_cex) / 1.45, n / 2 + 1), asp = 1, axes = FALSE, xlab = "", ylab = "")
  col <- pal[pmax(1, pmin(length(pal), findInterval(v, brk, all.inside = TRUE)))]
  a <- pos[Pairs$i]; b <- pos[Pairs$j]; x <- (a + b) / 2; y <- abs(b - a) / 2
  for (k in seq_along(v)) polygon(x[k] + c(-.5, 0, .5, 0), y[k] + c(0, .5, 0, -.5), col = col[k], border = NA)
  if (!is.null(mark)) points(x[mark], y[mark], pch = 16, cex = 0.12, col = "grey15")   # small dot = pairs to flag (FDR < 0.05)
  for (k in which(Pairs$i == kc | Pairs$j == kc)) polygon(x[k] + c(-.5, 0, .5, 0), y[k] + c(0, .5, 0, -.5), col = NA, border = "black", lwd = 0.6)
  text(pos[A$Shown], -0.6, A$Txt[A$Shown], srt = -55, adj = 0, cex = lab_cex, xpd = NA, col = ifelse(L$Loci[A$Shown] == 0, "red3", "black"))
  segments(pos[kc], -0.4, pos[kc], 0, col = "red3", lwd = 1.5)
  xs <- seq(n * 0.03, n * 0.33, length.out = length(pal) + 1)                 # key in the empty upper-left corner
  rect(head(xs, -1), n * 0.40, tail(xs, -1), n * 0.43, col = pal, border = NA)
  text(c(xs[1], mean(xs), tail(xs, 1)), n * 0.385, format(c(brk[1], 0, tail(brk, 1)), digits = 2), cex = key_cex, adj = c(0.5, 1))
  text(mean(xs), n * 0.45, key_lab, cex = key_cex * 1.1, adj = c(0.5, 0))
}
excess_palette <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))(21)

### Write a figure as PDF (vector) and 600-dpi PNG
save_fig <- function(file, draw, w = 7, h = 5) {
  pdf(paste0(file, ".pdf"), width = w, height = h, useDingbats = FALSE); draw(); dev.off()
  for (k in 1:5) {                         # a file briefly locked by a viewer or sync client: wait and retry
    ok <- tryCatch({ suppressWarnings(png(paste0(file, ".png"), width = w, height = h, units = "in", res = 600)); TRUE }, error = function(e) FALSE)
    if (ok) break; Sys.sleep(2)
  }
  if (!ok) stop("cannot write ", file, ".png (file locked?)")
  draw(); dev.off()
}

