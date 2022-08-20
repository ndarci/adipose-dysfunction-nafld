## collectMetadata.R
## consolidate metadata for all my DE genes into one table

library(caroline)
library(ggplot2)
library(ggthemes)
library(plyr)

# setwd('~/src/adiDys/bin/scripts')

# read in my DE genes
fibgenes <- read.table("../../data/DEgenes_fibrosis_adipose.txt")
stegenes <- read.table("../../data/DEgenes_steatosis_adipose.txt")
diagenes <- read.table("../../data/DEgenes_diagnosis_adipose.txt")
# infgenes <- read.table("../../data/DEgenes_inflammation_adipose.txt")
# balgenes <- read.table("../../data/DEgenes_ballooning_adipose.txt")
# unhgenes <- read.table("../../data/DEgenes_unhealthy_adipose.txt")
# stegenes_m <- read.table("../../data/DEgenes_steatosis_male.txt")
# stegenes_fm <- read.table("../../data/DEgenes_steatosis_female.txt")
# sexgenes <- read.table("../../data/DEgenes_sex.txt")
# stesexgenes <- read.table("../../data/DEgenes_steatosis_sex_interaction.txt")

# add group columns
fibgenes["group"] <- "fibrosis"
stegenes["group"] <- "steatosis"
diagenes["group"] <- "diagnosis"
# infgenes["group"] <- "inflammation"
# balgenes["group"] <- "ballooning"
# unhgenes["group"] <- "unhealthy"
# stegenes_m["group"] <- "steatosis_male"
# stegenes_fm["group"] <- "steatosis_female"
# sexgenes["group"] <- "sex"
# stesexgenes["group"] <- "steatosis_sex_interaction"

# combine into one dataframe
de_raw_sub <- rbind(fibgenes, stegenes, diagenes)#, infgenes, balgenes, unhgenes)
# de_raw <- rbind(fibgenes, stegenes, diagenes, infgenes, balgenes, unhgenes, stegenes_m, stegenes_fm, sexgenes, stesexgenes)
# de_raw <- rbind(fibgenes, stegenes, diagenes, infgenes, balgenes, unhgenes, stegenes_fm, sexgenes, stesexgenes)
# colnames(de_raw)[1] <- "Ensembl"
colnames(de_raw_sub)[1] <- "Ensembl"

# get all unique DE genes
uniqDEall <- unique(de_raw_sub$Ensembl)
write.table(uniqDEall, file = "../../data/allDEgenes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# find the number of total and unique DE genes for each group
grps <- unique(de_raw_sub$group)
de_grpByGene <- groupBy(de_raw_sub, by = "Ensembl", aggregation = c("max", "paste"), clmns = c("Ensembl", "group"))
de_uniqOnly <- de_grpByGene[-grep(",", de_grpByGene$group), ]
tots <- c()
uniqs <- c()
for (g in grps)
{
  totDE <- dim(de_raw_sub[de_raw_sub$group == g, ])[1]
  uniqDE <- dim(de_uniqOnly[de_uniqOnly$group == g, ])[1]
  tots <- c(tots, totDE)
  uniqs <- c(uniqs, uniqDE)
  fn <- paste("../../data/uniqDEgenes_", g, ".txt", sep = "")
  write.table(de_uniqOnly[de_uniqOnly$group == g, ]$Ensembl, file = fn, row.names = FALSE, col.names = FALSE, quote = FALSE)
}
deCounts <- data.frame(group = grps, total_DE_genes = tots, unique_DE_genes = uniqs)
deCounts["unique_pct"] <- deCounts$unique_DE_genes / deCounts$total_DE_genes
# plot unique percentage
ord <- factor(deCounts$group, levels = c("steatosis", "fibrosis", "diagnosis"))#"inflammation", "unhealthy", "ballooning"))
ggplot(data = deCounts, aes(x = ord, y = unique_pct)) + geom_bar(stat = "identity") + ylim(0,1) + xlab("group") + ggtitle("Percent unique DE genes")

# add secreted/not secreted column
# read in list of secreted proteins
secreted <- read.table("../../data/protein_class_Secreted_MDSEC.tsv", header = TRUE, sep = "\t", fill = TRUE)
# find my DE genes that match to secreted list, mark them as secreted
de_secr <- join(de_raw_sub, secreted, by = "Ensembl", type = "inner")
de_secr <- de_secr[ , c("Ensembl", "group")]
de_secr["secreted"] <- TRUE
# mark all genes that don't match to the list as NOT secreted
de <- de_raw_sub[ ! de_raw_sub$Ensembl %in% de_secr$Ensembl, ]
de["secreted"] <- FALSE
# put the secreted and not secreted back together into de table
de <- rbind(de, de_secr)
colnames(de)[1] <- "gene_ID"

# add PANTHER protein class column
# read in PANTHER data
panth <- read.table("../../data/pantherGeneList_all.txt", header = TRUE, sep = '\t', quote = "", na.strings = c("", "NA"))
# select desired PANTHER columns
panth <- panth[ , c("Mapped.IDs", "PANTHER.Protein.Class", "PANTHER.GO.Slim.Molecular.Function", 
                    "PANTHER.GO.Slim.Biological.Process", "PANTHER.GO.Slim.Cellular.Component",
                    "GO.database.MF.Complete", "GO.database.BP.Complete", "GO.database.CC.Complete")]
colnames(panth) <- c("gene_ID", "PANTHER_prot_class", "goMFslim", "goBPslim", "goCCslim", "goMF", "goBP", "goCC")
# merge onto de table
de <- join(de, panth, by = "gene_ID", type = "left")

# add column for UP/DOWN expression in sick tissue
# UP means higher expression in sick patients
# DOWN means higher expression in healthy patients
# import toptables from each DE test
TTfib <- read.table("../../data/topTable_fibrosis_adipose.txt", row.names = 1)
TTste <- read.table("../../data/topTable_steatosis_adipose.txt", row.names = 1)
TTdia <- read.table("../../data/topTable_diagnosis_adipose.txt", row.names = 1)
# TTinf <- read.table("../../data/topTable_inflammation_adipose.txt", row.names = 1)
# TTbal <- read.table("../../data/topTable_ballooning_adipose.txt", row.names = 1)
# TTunh <- read.table("../../data/topTable_unhealthy_adipose.txt", row.names = 1)
# TTste_m <- read.table("../../data/topTable_steatosis_male.txt", row.names = 1)
# TTste_fm <- read.table("../../data/topTable_steatosis_female.txt", row.names = 1)
# TTsex <- read.table("../../data/topTable_sex.txt", row.names = 1)
# TTsteSex <- read.table("../../data/topTable_steatosis_sex_interaction.txt", row.names = 1)
# add group column to keep track of where DE results came from
TTfib["group"] <- "fibrosis"
TTste["group"] <- "steatosis"
TTdia["group"] <- "diagnosis"
# TTinf["group"] <- "inflammation"
# TTbal["group"] <- "ballooning"
# TTunh["group"] <- "unhealthy"
# TTste_m["group"] <- "steatosis_male"
# TTste_fm["group"] <- "steatosis_female"
# TTsex["group"] <- "sex"
# TTsteSex["group"] <- "steatosis_sex_interaction"
# add gene_ID column for later merging
TTfib["gene_ID"] <- rownames(TTfib)
TTste["gene_ID"] <- rownames(TTste)
TTdia["gene_ID"] <- rownames(TTdia)
# TTinf["gene_ID"] <- rownames(TTinf)
# TTbal["gene_ID"] <- rownames(TTbal)
# TTunh["gene_ID"] <- rownames(TTunh)
# TTste_m["gene_ID"] <- rownames(TTste_m)
# TTste_fm["gene_ID"] <- rownames(TTste_fm)
# TTsex["gene_ID"] <- rownames(TTsex)
# TTsteSex["gene_ID"] <- rownames(TTsteSex)
# combine all the toptables
# allTT <- rbind(TTfib, TTste, TTdia, TTinf, TTbal, TTunh, TTste_m, TTste_fm, TTsex, TTsteSex)
# allTT <- rbind(TTfib, TTste, TTdia, TTinf, TTbal, TTunh, TTste_fm, TTsex, TTsteSex)
allTT <- rbind(TTfib, TTste, TTdia)#, TTinf, TTbal, TTunh)
# define new column and rename values to UP and DOWN
allTT["up_in_sick"] <- allTT$logFC > 0
allTT[allTT$up_in_sick == TRUE, "up_in_sick"] <- "UP"
allTT[allTT$up_in_sick == FALSE, "up_in_sick"] <- "DOWN"
# allTT <- allTT[, c("gene_ID", "up_in_sick", "logFC", "B", "adj.P.Val")]
# colnames(allTT) <- c("gene_ID", "up_in_sick", "logFC", "B", "DE_adjP")
colnames(allTT)[4] <- "DE_P"
colnames(allTT)[5] <- "DE_adjP"
# merge onto de table
de <- join(de, allTT, by = c("gene_ID", "group"), type = "left")

# read in filtered annotation file from GENCODE
annot <- read.table("../../data/gencodeV19_annotations_filtered.txt", header = TRUE, sep = '\t')
colnames(annot) <- c("gene_ID", "gene_ID_GENCODEver", "chromosome", "feature", "startPos", "endPos", "strand", "gene_type", "gene_symbol")
annot <- annot[ , c("gene_ID", "gene_ID_GENCODEver", "chromosome", "startPos", "endPos", "strand", "gene_type")]
de <- join(de, annot, by = "gene_ID", type = "left")

# add GTEx expression information
# read in GTEx expression data: 54 total tissue types
gtex <- read.table("../../data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", header = TRUE, skip = 2, sep = '\t')
colnames(gtex)[1] <- "gene_ID_GTExVer"
colnames(gtex)[2] <- "gene_symbol"
novid <- strsplit(gtex$gene_ID_GTExVer, ".", fixed = TRUE)
gtex["gene_ID"] <- sapply(novid, "[[", 1)
# # extract just the gene names from GTEx table
# gtex_names <- gtex[ , c("gene_ID", "gene_ID_GTExVer", "gene_symbol")]
# merge the names onto the de table
de <- join(de, gtex, by = "gene_ID", type = "left")
# add adipose to liver expression ratio column
de["adipose_over_liver_expr"] <- de$Adipose...Subcutaneous / de$Liver

# extract the GTEx expression information for just our DE genes
gtexDE <- join(de[c("gene_ID")], gtex, by = "gene_ID", type = "left")
gtexDE <- unique(na.omit(gtexDE))
# gc <- gc[-grep(",", gc$cellType_cluster), ]
gtexDE <- gtexDE[-grep("PAR_Y", gtexDE$gene_ID_GTExVer), ]
# rownames(gtexDE) <- gtexDE$gene_symbol
gtexDE <- subset(gtexDE, select = -c(gene_ID, gene_ID_GTExVer, gene_symbol))
write.table(gtexDE, file = "../../data/deGenesMedianTPM_GTEx.txt", sep = '\t', quote = FALSE)

# add DE in liver column
# read in liver DE data for same 6 tests as were done in adipose
livFib <- read.table("../../data_liver/topTable_fibrosis_liver.txt")
livSte <- read.table("../../data_liver/topTable_steatosis_liver.txt")
# livInf <- read.table("../../data_liver/topTable_inflammation_liver.txt")
livDia <- read.table("../../data_liver/topTable_diagnosis_liver.txt")
# livBal <- read.table("../../data_liver/topTable_ballooning_liver.txt")
# livUnh <- read.table("../../data_liver/topTable_unhealthy_liver.txt")
# add group and ID columns
livFib["liver_group"] <- "fibrosis"
livFib["gene_ID_ver"] <- rownames(livFib)
livSte["liver_group"] <- "steatosis"
livSte["gene_ID_ver"] <- rownames(livSte)
# livInf["liver_group"] <- "inflammation"
# livInf["gene_ID_ver"] <- rownames(livInf)
livDia["liver_group"] <- "diagnosis"
livDia["gene_ID_ver"] <- rownames(livDia)
# livBal["liver_group"] <- "ballooning"
# livBal["gene_ID_ver"] <- rownames(livBal)
# livUnh["liver_group"] <- "unhealthy"
# livUnh["gene_ID_ver"] <- rownames(livUnh)
liv <- rbind(livFib, livSte, livDia)#, livInf, livBal, livUnh)
# filter for significant DE genes
liv <- liv[liv$adj.P.Val < 0.05, ]
# remove version numbers
splLiv <- strsplit(liv$gene_ID_ver, ".", fixed = TRUE)
liv["gene_ID"] <- sapply(splLiv, "[[", 1)
# save table of all liver results together
write.table(liv, "../../data_liver/deGenes_all_liver.txt", quote = F, row.names = F, sep = '\t')
# group on gene ID
liv <- liv[ , c("gene_ID", "liver_group")]
liv <- groupBy(liv, by = "gene_ID", aggregation = c("max", "paste"), clmns = c("gene_ID", "liver_group"))
# add columns to the adipose table
de["DE_in_liver"] <- de$gene_ID %in% liv$gene_ID
de <- join(de, liv, by = "gene_ID", type = "left")

# add cell type marker information
# import data
ctm <- read.table("../../data/SC-Sub_Minna_snRNA20200617_DIEM.demux.noMT.pc_markers_CellTySingleR.Clust_fdr0.05.txt", sep = '\t', header = TRUE)
# filter for significance
ctm <- ctm[ctm$p_val_adj < 0.05, ]
# select the columns we want and rename them
gc <- ctm[ , c("gene", "marker.clust", "p_val_adj")]
colnames(gc) <- c("gene_symbol", "cellType_cluster", "cellType_adjP")
# # select only unique cell type marker genes
gc <- groupBy(gc, by = "gene_symbol", aggregation = c("max", "paste", "mean"), clmns = c("gene_symbol", "cellType_cluster", "cellType_adjP"))
# gc <- gc[-grep(",", gc$cellType_cluster), ]
# merge onto de table
de <- join(de, gc, by = "gene_symbol", type = "left")
# count the number of unique cell type markers out of all my DE genes
de_ctm <- unique(de[, c("gene_ID", "cellType_cluster")])
# how many of my genes are unique cell type markers?
tctm <- table(is.na(de_ctm$cellType_cluster))
ctmPct <- tctm[[1]] / (tctm[[1]] + tctm[[2]])
#ctmPct
# how many genes in each cluster?
# table(de_ctm$cellType_cluster)

# # grab some novel transcripts
# novTrans <- de[is.na(de$gene_symbol), ]
# novTransIDs <- unique(novTrans$gene_ID)
# write.table(novTransIDs, file = "../../data/novelDEtranscripts.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# what columns do we want to focus on?
desfacs <- c("group", "gene_ID", "gene_symbol", "secreted", "up_in_sick", "DE_in_liver", 
             "Adipose...Subcutaneous", "Liver", "adipose_over_liver_expr", 
             "logFC", "B", "t", "AveExpr", "DE_P", "DE_adjP")
deSub <- de[ , desfacs]
# drop duplicate rows
de <- unique(de)
deSub <- unique(deSub)
deSub <- deSub[order(-deSub$adipose_over_liver_expr), ]

# save table of all DE genes
write.table(deSub, "../../data/metTable_all_subset.txt", sep = '\t', quote = F, row.names = F)

# highlight serum biomarker candidates
TPMthresh <- 30
ratioThresh <- 10
deSub_sbc <- deSub[deSub$group %in% c("fibrosis", "steatosis", "diagnosis"),]#, "inflammation", "ballooning", 
                                      # "unhealthy"), ]
deSub_sbc = na.omit(deSub_sbc[,desfacs])

# check how filters would affect DE genes I want for PRS
for (prs_pheno in c("steatosis", "fibrosis", "diagnosis")) { #"inflammation", "ballooning", "unhealthy")) {
	df = deSub_sbc[deSub_sbc$group == prs_pheno,]
	print(paste0(prs_pheno, ": ", length(unique(df$gene_ID)), " DE genes"))
	print("NOT DE in liver")
	print(table(!(unique(df[,c("gene_ID", "DE_in_liver")])$DE_in_liver)))
	print("Adipose expression over threshold")
	print(table(unique(df[,c("gene_ID", "Adipose...Subcutaneous")])$Adipose...Subcutaneous > TPMthresh))
	print("Adipose:liver expression over threshold")
	print(table(unique(df[,c("gene_ID", "adipose_over_liver_expr")])$adipose_over_liver_expr > ratioThresh))
	print("Secreted")
	print(table(unique(df[,c("gene_ID", "secreted")])$secreted))
	df_filt = df[df$DE_in_liver == FALSE,]
	nfilt = length(unique(df_filt$gene_ID))
	print(paste0("Number passing all filters I care about: ", nfilt))
	write.table(unique(df_filt$gene_ID), 
		paste0("../../data/DE_genes_filter_except_secreted_", prs_pheno, "_grp.txt"),
		row.names = F, col.names = F, quote = F)
}

print(paste0("Starting genes: ", length(unique(deSub_sbc$gene_ID))))
deSub_sbc <- deSub_sbc[deSub_sbc$DE_in_liver == FALSE, ] # needs to NOT be DE in liver (all are DE in adipose)
print(paste0("After DE in liver filter: ", length(unique(deSub_sbc$gene_ID))))
deSub_sbc <- deSub_sbc[deSub_sbc$`Adipose...Subcutaneous` > TPMthresh, ] # needs to have detectable TPM
print(paste0("After adipose expression filter: ", length(unique(deSub_sbc$gene_ID))))
deSub_sbc <- deSub_sbc[deSub_sbc$adipose_over_liver_expr > ratioThresh, ] # needs to have significantly high adipose:liver expression
print(paste0("After adipose:liver expression filter: ", length(unique(deSub_sbc$gene_ID))))
deSub_sbc <- deSub_sbc[deSub_sbc$secreted == TRUE, ] # needs to be secreted
print(paste0("After secreted filter (final SBC count): ", length(unique(deSub_sbc$gene_ID))))
# deSub_sbc <- groupBy(deSub_sbc, by = "gene_symbol", aggregation = c("paste", rep("max", length(desfacs) - 1)), clmns = colnames(deSub_sbc))
write.table(deSub_sbc, file = "../../data/metTable_subset_serumBiomarkers.txt", sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

# make a new pretty and simple SBC metTable
library(tidyverse)
sbc.clean = deSub_sbc
sbc.clean['group'] = recode(sbc.clean$group, steatosis = 'Steatosis',
												fibrosis = 'Fibrosis',
												diagnosis = 'NASH')
sbc.clean = sbc.clean %>% 
					group_by(gene_symbol) %>%
					summarise(de_group = paste(group, collapse = ', '),
										secreted = secreted[1],
										DE_in_liver = DE_in_liver[1],
										adipose_median_tpm = Adipose...Subcutaneous[1],
										adipose_over_liver_expr = adipose_over_liver_expr[1])
colnames(sbc.clean) = c('Gene', 'Adipose DE traits', 'Secreted', 'DE in liver',
												'Adipose expression (median TPM)', 'Adipose:liver expression')
write.table(sbc.clean, file = '../../data/metTable_subset_serumBiomarkers_clean.txt', row.names = F)

sbc <- unique(deSub_sbc$gene_symbol)
sbc_ensembl <- unique(deSub_sbc$gene_ID)
write.table(sbc, file = "../../data/sbc_IDs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sbc_ensembl, file = "../../data/sbc_IDs_ensembl.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# # save TMM-normalized CPMs (just SBCs)
# cpm <- read.csv("../../data/TMMnorm_cpms.csv")
# cpm_sbc <- cpm[cpm$X %in% deSub_sbc$gene_ID,]
# colnames(cpm_sbc)[1] <- "gene_ID"
# cpm_sbc <- join(cpm_sbc, deSub_sbc[,c("gene_ID", "gene_symbol")], by = "gene_ID", type = "inner")
# cpm_sbc <- unique(cpm_sbc)
# rownames(cpm_sbc) <- cpm_sbc$gene_symbol
# cpm_sbc <- subset(cpm_sbc, select = -c(gene_ID, gene_symbol))
# write.csv(cpm_sbc, "../../data/TMMnorm_cpms_SBC.csv", row.names = T)

# # plot adipose:liver expression ratios in SBCs
# theme_set(theme_bw())
# ord2 <- deSub_sbc[order(-deSub_sbc$adipose_over_liver_expr), ]$gene_symbol
# ratioBoxPlot <- ggplot(data = deSub_sbc, aes(x = factor(gene_symbol, level = ord2), y = adipose_over_liver_expr)) + 
#   geom_col() + 
#   geom_text(aes(label = floor(adipose_over_liver_expr)), position = position_dodge(width = 1), vjust = -0.25, size = 10) +
#   xlab("Gene") + ylab("Adipose over liver expression") +
#   theme(text = element_text(family = "Times New Roman"),
#         axis.text.x = element_text(size = 30, angle = 45, vjust = 1, hjust = 1),
#         axis.text.y = element_text(size = 30),
#         axis.title = element_text(size = 40),
#         axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
#         plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         panel.grid = element_blank())
# ratioBoxPlot
# ggsave(plot = ratioBoxPlot, filename = "figs/adiposeOverLiverBoxplot.png", width = 19.875, height = 13.5, units = "in")

# separate out by group, with the interesting columns we want
deFibSub <- deSub[deSub$group == "fibrosis", ]
deSteSub <- deSub[deSub$group == "steatosis", ]
deDiaSub <- deSub[deSub$group == "diagnosis", ]
# deInfSub <- deSub[deSub$group == "inflammation", ]
# deBalSub <- deSub[deSub$group == "ballooning", ]
# deUnhSub <- deSub[deSub$group == "unhealthy", ]
# deSteSub_m <- deSub[deSub$group == "steatosis_male", ]
# deSteSub_fm <- deSub[deSub$group == "steatosis_female", ]
# deSexSub <- deSub[deSub$group == "sex", ]
# deSteSexSub <- deSub[deSub$group == "steatosis_sex_interaction", ]

# add p-value and logFC
desTTfacs <- c("gene_ID", "logFC", "adj.P.Val")
deFibSub <- join(deFibSub, TTfib[ , desTTfacs], by = "gene_ID")
deSteSub <- join(deSteSub, TTste[ , desTTfacs], by = "gene_ID")
# deInfSub <- join(deInfSub, TTinf[ , desTTfacs], by = "gene_ID")
deDiaSub <- join(deDiaSub, TTdia[ , desTTfacs], by = "gene_ID")
# deBalSub <- join(deBalSub, TTbal[ , desTTfacs], by = "gene_ID")
# deUnhSub <- join(deUnhSub, TTunh[ , desTTfacs], by = "gene_ID")
# deSteSub_m <- join(deSteSub_m, TTste_m[ , desTTfacs], by = "gene_ID")
# deSteSub_fm <- join(deSteSub_fm, TTste_fm[ , desTTfacs], by = "gene_ID")
# deSexSub <- join(deSexSub, TTsex[ , desTTfacs], by = "gene_ID")
# deSteSexSub <- join(deSteSexSub, TTsteSex[ , desTTfacs], by = "gene_ID")

# # which serum biomarker candidates are sex-specific?
# deSexSub[deSexSub$gene_symbol %in% sbc, ]$gene_symbol
# deSteSub_fm[deSteSub_fm$gene_symbol %in% sbc, ]$gene_symbol
# deSteSexSub[deSteSexSub$gene_symbol %in% sbc, ]$gene_symbol

# write out metadata tables (interesting columns only)
write.table(deFibSub, file = "../../data/metTable_fibrosis_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(deSteSub, file = "../../data/metTable_steatosis_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(deDiaSub, file = "../../data/metTable_diagnosis_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(deInfSub, file = "../../data/metTable_inflammation_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(deBalSub, file = "../../data/metTable_ballooning_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(deUnhSub, file = "../../data/metTable_unhealthy_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(deSteSub_m, file = "../../data/metTable_steatosis_subset_male.txt", sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(deSexSub, file = "../../data/metTable_sex_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(deSteSexSub, file = "../../data/metTable_steatosis_sex_interaction_subset.txt", sep = '\t', quote = FALSE, row.names = FALSE)

