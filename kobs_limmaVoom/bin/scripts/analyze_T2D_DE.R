# analyze_T2D_DE.R
# read in the DE results for type 2 diabetes and check for the SBCs

setwd("~/src/adiDys/bin/scripts/")

# read in T2D toptable from limma-voom
tt <- read.table("../../data/topTable_T2D_adipose.txt")
tt_sig <- tt[tt$adj.P.Val < 0.05,]

# read in list of SBCs
sbc <- read.table("../../data/metTable_subset_serumBiomarkers.txt", sep = '\t', header = T)

# find SBCs in T2D results
tt_sbc <- tt[rownames(tt) %in% sbc$gene_ID,]
tt_sbc["gene_ID"] <- rownames(tt_sbc)
tt_sbc <- merge(tt_sbc, unique(sbc[,c("gene_ID", "gene_symbol")]), by = "gene_ID")
rownames(tt_sbc) <- tt_sbc$gene_symbol
tt_sbc <- tt_sbc[2:7]
