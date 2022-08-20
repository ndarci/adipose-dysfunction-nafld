library(tidyverse)
library(ggplot2)
library(ggpubr)

time_order = c("Baseline", "24h", "4D", "7D")
sbc_order = c("CCDC80", "SOD3")

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

read_in_de <- function(nominal = F) {
	# import DE results
	de = data.frame()
	for (tp in time_order) {
		for (sbc in sbc_order) {
			# for (covariate in c(paste0("SV1_", sbc), "unique_map_pct", "")) {
			# 	if (covariate == "") { 
					covcorrect = F 
					covariate = ""
				# } else { covcorrect = T }
				fn = paste0("../../data/topTable_", 
							tp, "_", 
							sbc, 
							rep(paste0("_covcorrected_", covariate), covcorrect), 
							rep("_nominal", nominal),
							".txt")
				newdata = read.table(fn, header = T, sep = '\t')
				if (nrow(newdata) > 0) {
					newdata["timepoint"] = tp
					newdata["knockdown"] = sbc
					newdata["covcorrect"] = covcorrect
					newdata["covariate"] = covariate
					de = rbind(de, newdata)
				# }
			}
		}
	}
	return(de)
}

# read in all the de results... ignore covcorrect analyses
de = read_in_de()
de['gene_ID_ver'] = de$gene_ID
de['gene_ID'] = stripVersion(de$gene_ID_ver)

# define up in knockdown column
de['up_in_KD'] = de$logFC > 0

# count DE genes per condition
de_count = de %>% 
	group_by(timepoint, knockdown) %>%
	summarise(de_genes_SV1 = length(gene_ID[grepl("SV1", covariate)]),
				de_genes_nocov = length(gene_ID[covariate == ""]),
				de_genes_ump = length(gene_ID[covariate == "unique_map_pct"]),
				de_genes_shared = length(intersect(intersect(gene_ID[covariate == ""], 
										gene_ID[grepl("SV1", covariate)]),
										gene_ID[covariate == "unique_map_pct"])))
# print(de_count)
write.table(de_count, "../../data/numDEgenes_allconditions.txt", quote = F, row.names = F, sep = '\t')

# plot DE genes per condition, and overlap across covcorrect/no covcorrect
for (sbc in c("CCDC80", "SOD3")) {
	de_count_sub = de_count[de_count$knockdown == sbc,] %>%
				pivot_longer(cols = starts_with("de_genes"),
							names_to = "de_type",
							names_prefix = "de_genes_",
							values_to = "num_de")
	de_count_sub$de_type = factor(de_count_sub$de_type, level = c("nocov", "SV1", "ump", "shared"))
	de_count_sub$timepoint = factor(de_count_sub$timepoint, level = c("Baseline", "24h", "4D", "7D"))
	numdeplot = ggbarplot(de_count_sub, x = "timepoint", y = "num_de",
				fill = "de_type", position = position_dodge(0.8)) +
			xlab("Timepoint") + 
			ylab(paste0("Number of DE genes in ", sbc, " knockdown")) +
			theme(legend.position = "right")  
			# scale_fill_discrete(name = "", labels = c())
	ggsave(paste0("../../figs/numDEgenes_", sbc, "_allconditions_bar.png"), numdeplot)
}

# for (sbc in sbc_order) {
# 	for (c in c(T, F)) {
# 		de_count_sub = de_count[de_count$knockdown == sbc &
# 						de_count$covcorrect == c,] %>%
# 				arrange(factor(timepoint, levels = time_order))
# 		fn = paste0("../../data/numDEgenes_", sbc, rep("_covcorrected", c), ".txt")
# 		write.table(de_count_sub, fn, quote = F, row.names = F, sep = '\t')
# 	print(de_count_sub)
# 	}
# }

# search for all SBCs in these results
# sbc = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/sbc_IDs.txt")$V1

# de_sbc = de[de$gene_name %in% sbc,]
# write.table(de_sbc[de_sbc$covcorrect == T,], "../../data/knockdown_de_overlap_SBCs_covcorrected.txt", quote = F, sep = '\t', row.names = F)
# write.table(de_sbc[de_sbc$covcorrect == F,], "../../data/knockdown_de_overlap_SBCs.txt", quote = F, sep = '\t', row.names = F)

# # search for PPARG and ADIPOQ in the results
# de_adip = de[de$gene_name %in% c("PPARG", "ADIPOQ"),]
# write.table(de_adip, "../../data/knockdown_de_overlap_adipogenes.txt", quote = F, sep = '\t', row.names = F)





# verify knockdown with nominal p-values
# read in data
de_nom = read_in_de(nominal = T)
de_nom_kd = de_nom[de_nom$gene_name %in% sbc_order,]

# filter for significant NOMINAL p 
de_nom_kd = de_nom_kd[de_nom_kd$P.Value < 0.05,]

# take only the same gene that was knocked down
de_nom_kd = de_nom_kd[de_nom_kd$gene_name == de_nom_kd$knockdown,]
write.table(de_nom_kd, "../../data/knockdown_verification_nomPval.txt", row.name = F, quote = F)







# search for adipogenesis genes in the results
adigen = read.table("../../data/adipogenesis_genelist_ensembl.txt")$V1
adipocyte = read.table("../../data/adipocyte_genelist_ensembl.txt")$V1
preadipocyte = read.table("../../data/preadipocyte_genelist_ensembl.txt")$V1
adigen_adi_preadi = c(adigen, adipocyte, preadipocyte)

# match genes with adipogenesis pathway annotations
inhibdiff = data.frame(gene_name = c("DLK1", "GATA2", "GATA3", "GATA4", "WNT1", "WNT5B", "WNT10B", "NDN"),
						wiki_annot = "Inhibitors of preadipoctye to adipoctye diff.")
tf = data.frame(gene_name = c("CEBPA", "CEBPB", "CEBPD", "PPARA", "PPARG", "PPARD", "RXRA", "RXRG", "RARA",
		"NR1H3", "CREB1", "MEF2A", "MEF2B", "MEF2C", "MEF2D", "SREBF1", "MBNL1", "NCOA1", 
		"NR2F1", "RORA", "NRIP1", "NCOR1", "NCOR2"), 
		wiki_annot = "TFs/modulators")
growth = data.frame(gene_name = c("INS", "IGF1", "NR3C1", "GH1", "TGFB1", "PRLR"),
		wiki_annot = "Growth factors/hormones")
fulldiff = data.frame(gene_name = c("GTF3A", "FAS", "SLC2A4", "PLIN", "LPL", "PPARGC1A", "UCP1", "LIPE"),
		wiki_annot = "Fully diff. adipocyte markers")
insul = data.frame(gene_name = c("IRS1", "IRS2", "IRS3P", "IRS4"),
		wiki_annot = "Insulin action genes")
secrete = data.frame(gene_name = c("LEP", "ADIPOQ", "IL6", "TNF", "RETN", "DF", "ADFP", "ADPN", "AGT", "PBEF1", "SPOCK"),
		wiki_annot = "Adipocyte secretory products")
lipodys = data.frame(gene_name = c("BSCL2", "LMNA", "AGPAT2", "ZMPSTE24", "LPIN1", "LPIN2", "LPIN3"),
		wiki_annot = "Possible lipodystrophy genes")
adigen_df = rbind(inhibdiff, tf, growth, fulldiff, insul, secrete, lipodys)
othergenes = unique(de[de$gene_ID %in% adigen &
				!(de$gene_name %in% adigen_df$gene_name),]$gene_name)
other = data.frame(gene_name = othergenes,
		wiki_annot = "Miscellaneous factors")
adigen_df = rbind(adigen_df, other)
adigen_df["adipogenesis"] = T

de = de[order(de$adj.P.Val, decreasing = F),]

# merge on adipogenesis annotations
de = merge(de, adigen_df, by = "gene_name", all.x = T)
de[is.na(de$adipogenesis),"adipogenesis"] <- F

# save genes in adipogenesis wikipathway for sod3 and ccdc80
for (sbc in sbc_order) {
	de_wiki = de[de$knockdown == sbc & de$adipogenesis == T,]
	write.table(de_wiki, paste0('../../data/wikipathway_overlap_', sbc, '.txt'), row.names = F)
}

# merge on adipocyte and preadipocyte information
de["adipocyte"] = de$gene_ID %in% adipocyte
de["preadipocyte"] = de$gene_ID %in% preadipocyte

# add secreted/not secreted column
# read in list of secreted proteins
secreted <- read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/protein_class_Secreted_MDSEC.tsv", header = TRUE, sep = "\t", fill = TRUE)
# find my DE genes that match to secreted list, mark them as secreted
de['secreted'] = de$gene_ID %in% secreted$Ensembl

# add PANTHER protein class column
# read in PANTHER data
panth <- read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/pantherGeneList_all.txt", header = TRUE, sep = '\t', quote = "", na.strings = c("", "NA"))
# select desired PANTHER columns
panth <- panth[ , c("Mapped.IDs", "PANTHER.Protein.Class", "PANTHER.GO.Slim.Molecular.Function", 
                    "PANTHER.GO.Slim.Biological.Process", "PANTHER.GO.Slim.Cellular.Component")]
colnames(panth) <- c("gene_ID", "PANTHER_prot_class", "goMFslim", "goBPslim", "goCCslim")
# merge onto de table
de <- merge(de, panth, by = "gene_ID", all.x = T)


# add gtex data for adipose and liver
gtex <- read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", header = TRUE, skip = 2, sep = '\t')
colnames(gtex)[1] <- "gene_ID_GTExVer"
colnames(gtex)[2] <- "gene_symbol"
gtex['gene_ID'] = stripVersion(gtex$gene_ID_GTExVer)
gtex = gtex[,c('gene_ID', 'Adipose...Subcutaneous', 'Adipose...Visceral..Omentum.', 'Liver')]
# merge the names onto the de table
de <- merge(de, gtex, by = 'gene_ID', all.x = T)
# add adipose to liver expression ratio column
de["adipose_over_liver_expr"] <- de$Adipose...Subcutaneous / de$Liver

# look at each important annotation specifically
annot_map = list(`TFs/modulators` = "TF_modulator",
				`Fully diff. adipocyte markers` = "fully_diff_adipocyte", 
				`Adipocyte secretory products` = "adipocyte_secrete", 
				`Inhibitors of preadipoctye to adipoctye diff.` = "inhib_adipogen")
for (annot in c("TFs/modulators", 
				"Fully diff. adipocyte markers", 
				"Adipocyte secretory products", 
				"Inhibitors of preadipoctye to adipoctye diff.")) {
	de_adi_annot = de[de$wiki_annot == annot & !is.na(de$wiki_annot),]
	write.table(de_adi_annot, paste0("../../data/metTable_allgroups_annot_", annot_map[annot], ".txt"), row.names = F, quote = F, sep = '\t')
}

# find which DE genes show up in the srebf1 wikipathway list
# read in the pathway gene list
path = read.table('../../data/srebf1_wikipathway_ensembl.txt')$x
de['SREBF1_wikipathway'] = de$gene_ID %in% path

# find which show up in the mayaan lab srebf1 target list
hzm = read.table('../../data/srebf1_targetgenes_mayaanlab.txt')$V1 
de['SREBF1_mayaan'] = de$gene_name %in% hzm 

# check where WGCNA module TFs show up in knockdown
aly.tf = c('ELF4', 'GTF2E2', 'TCF23', 'IRX6', 'ZBTB7B', 'REPIN1')
lsb.tf = c('ID2', 'ETV3', 'ZNF281', 'TEAD1', 'ELOA')
de['A_lightyellow_tf'] = de$gene_name %in% aly.tf
de['L_saddlebrown_tf'] = de$gene_name %in% lsb.tf

# write out results by DE test
de_out = de[,!(colnames(de) %in% c("covcorrect", "covariate", "gene_ID_ver"))]
for (sbc in sbc_order) {
	write.table(de_out[(de_out$SREBF1_wikipathway == T |
							de_out$SREBF1_mayaan == T) &
						de_out$knockdown == sbc,],
				paste0("../../data/metTable_srebf1_only_", sbc, ".txt"))
	for (tp in time_order) {
		de_out_sub = de_out[de_out$timepoint == tp &
							de_out$knockdown == sbc,]
		fn = paste0("../../data/metTable_", tp, "_", sbc, ".txt")
		write.table(de_out_sub, fn, row.names = F, quote = F, sep = '\t')
	}
}
