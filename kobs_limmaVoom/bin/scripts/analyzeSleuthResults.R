## analyzeSleuthResults.R
# import sleuth results from hoffman and compare them to limma results

library(plyr)

setwd("~/src/adiDys/bin/scripts")

# import transcript DE results
# tx <- read.table(file = "../../data/sleuthResultsTranscripts.txt", header = TRUE)
tx <- read.table(file = "/u/home/n/nikodm/project-pajukant/kobs_kallisto/data/sleuthResultsTranscripts.txt", header = T)

# extract versioned/non-versioned ENST and ENSG IDs
splitme <- strsplit(tx$target_id, "|", fixed = TRUE)
enst_ver <- sapply(splitme, "[[", 1)
ensg_ver <- sapply(splitme, "[[", 2)
enst <- sapply(strsplit(enst_ver, ".", fixed = TRUE), "[[", 1)
ensg <- sapply(strsplit(ensg_ver, ".", fixed = TRUE), "[[", 1)
tx["transcript_id_ver"] <- enst_ver
tx["transcript_id"] <- enst
tx["gene_id_ver"] <- ensg_ver
tx["gene_id"] <- ensg

# import annotations to get gene symbol <-> ID mappings
annot <- read.table("../../data/gencodeV26_annotations_filtered.txt", header = TRUE, sep = '\t')
colnames(annot) <- c("gene_id", "gene_id_GENCODEver", "chromosome", "feature", "startPos", "endPos", "strand", "gene_type", "gene_symbol")
annot <- annot[ , c("gene_id", "gene_id_GENCODEver", "chromosome", "startPos", "endPos", "strand", "gene_type", "gene_symbol")]

# merge annotations onto DE results
tx <- join(tx, annot, by = "gene_id", type = "left")

# check which SBCs from the limma results are DE according to sleuth
sbc <- read.table("../../data/metTable_subset_serumBiomarkers.txt", header = TRUE, sep = "\t")
sbcIDs <- sbc$gene_symbol
tx_sbc <- tx[tx$gene_symbol %in% sbcIDs, ]

# get significant DE transcripts
tx_sig <- tx[tx$qval < 0.05,]

# check which of these transcripts come from genes with multiple isoforms
# output list of genes that DE transcripts come from
write.table(unique(tx_sig$gene_id), file = "../../data/DEtranscriptGenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
# plug that^ gene list into ENSEMBL GRCh37 BioMart database, and have it return Gene ID/TX ID
bm <- read.table("../../data/DEtranscriptGenes_withTX_BioMart.txt", header = TRUE, sep = "\t")
# count each gene's number of isoforms, and merge it back onto the transcripts
numtx <- data.frame(table(bm$Gene.stable.ID))
colnames(numtx) <- c("gene_id", "num_isoforms")
tx_sig <- join(tx_sig, numtx, by = "gene_id", type = "left")

# what proportion of DE transcripts come from genes with multiple isoforms?
mean(tx_sig$num_isoforms > 1)






