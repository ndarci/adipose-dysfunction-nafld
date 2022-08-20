# look for the genes/transcripts assoc with isoQTLs in my DE results

# import all adipose DE
TTfib <- read.table("../../data/topTable_fibrosis_adipose.txt", row.names = 1)
TTste <- read.table("../../data/topTable_steatosis_adipose.txt", row.names = 1)
TTdia <- read.table("../../data/topTable_diagnosis_adipose.txt", row.names = 1)
TTinf <- read.table("../../data/topTable_inflammation_adipose.txt", row.names = 1)
TTbal <- read.table("../../data/topTable_ballooning_adipose.txt", row.names = 1)
TTunh <- read.table("../../data/topTable_unhealthy_adipose.txt", row.names = 1)
# add group column to keep track of where DE results came from
TTfib["group"] <- "fibrosis"
TTste["group"] <- "steatosis"
TTdia["group"] <- "diagnosis"
TTinf["group"] <- "inflammation"
TTbal["group"] <- "ballooning"
TTunh["group"] <- "unhealthy"
TTfib["gene_ID"] <- rownames(TTfib)
TTste["gene_ID"] <- rownames(TTste)
TTdia["gene_ID"] <- rownames(TTdia)
TTinf["gene_ID"] <- rownames(TTinf)
TTbal["gene_ID"] <- rownames(TTbal)
TTunh["gene_ID"] <- rownames(TTunh)
adi <- rbind(TTfib, TTste, TTdia, TTinf, TTbal, TTunh)
adi = subset(adi, adi$adj.P.Val < 0.05)

# import all 6 liver DE tables
livFib <- read.table("../../data_liver/topTable_fibrosis_liver.txt")
livSte <- read.table("../../data_liver/topTable_steatosis_liver.txt")
livInf <- read.table("../../data_liver/topTable_inflammation_liver.txt")
livDia <- read.table("../../data_liver/topTable_diagnosis_liver.txt")
livBal <- read.table("../../data_liver/topTable_ballooning_liver.txt")
livUnh <- read.table("../../data_liver/topTable_unhealthy_liver.txt")
# add group and ID columns
livFib["liver_group"] <- "fibrosis"
livFib["gene_ID_ver"] <- rownames(livFib)
livSte["liver_group"] <- "steatosis"
livSte["gene_ID_ver"] <- rownames(livSte)
livInf["liver_group"] <- "inflammation"
livInf["gene_ID_ver"] <- rownames(livInf)
livDia["liver_group"] <- "diagnosis"
livDia["gene_ID_ver"] <- rownames(livDia)
livBal["liver_group"] <- "ballooning"
livBal["gene_ID_ver"] <- rownames(livBal)
livUnh["liver_group"] <- "unhealthy"
livUnh["gene_ID_ver"] <- rownames(livUnh)
liv <- rbind(livFib, livSte, livInf, livDia, livBal, livUnh)
# filter for significant DE genes
liv <- liv[liv$adj.P.Val < 0.05, ]
# remove version numbers
splLiv <- strsplit(liv$gene_ID_ver, ".", fixed = TRUE)
liv["gene_ID"] <- sapply(splLiv, "[[", 1)

# import huiling's isoQTL genes
# note: had to add a newline at the very end of this file
iso = read.table("../../data/isoQTL_transcripts.txt", sep = "|")
colnames(iso) = c("tx_ID_ver", "gene_ID_ver", "gene_sym_ver", "gene_sym", "x", "gene_type", "x2")
# remove version numbers
splIso <- strsplit(iso$gene_ID_ver, ".", fixed = TRUE)
iso["gene_ID"] <- sapply(splIso, "[[", 1)
splIso2 <- strsplit(iso$tx_ID_ver, ".", fixed = TRUE)
iso["tx_ID"] <- sapply(splIso2, "[[", 1)

# check for huiling's DE genes in my results
ovLiv = merge(liv, iso, by = "gene_ID", all.x = F, all.y = F)
ovAdi = merge(adi, iso, by = "gene_ID", all.x = F, all.y = F)

write.table(ovLiv, "../../data/isoQTL_DE_gene_overlap_liver.txt", quote = F, row.names = F)
write.table(ovAdi, "../../data/isoQTL_DE_gene_overlap_adipose.txt", quote = F, row.names = F)

# import tx-level DE results
slu = read.table("/u/project/pajukant/nikodm/kobs_kallisto/data/sleuthResultsTranscripts.txt", header = T)
slu = subset(slu, slu$qval < 0.05)

# get gene and tx IDs
spltar = strsplit(slu$target_id, "|", fixed = T)
slu["tx_ID_ver"] = sapply(spltar, "[[", 1)
slu["gene_ID_ver"] = sapply(spltar, "[[", 2)
stx = strsplit(slu$tx_ID_ver, ".", fixed = T)
slu["tx_ID"] = sapply(stx, "[[", 1)
sg = strsplit(slu$gene_ID_ver, ".", fixed = T)
slu["gene_ID"] = sapply(sg, "[[", 1)

# check for DE tx's in tx results
ovSlu = merge(slu, iso, by = "tx_ID")









