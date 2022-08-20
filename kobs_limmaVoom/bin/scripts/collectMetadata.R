## collectMetadata.R
## consolidate metadata for all my DE genes into one table
library(tidyverse)
# library(venn)
# library(ggvenn)
# library(gridExtra)

hist_order = c('steatosis', 'fibrosis', 'diagnosis')

# function to read in kobs DE results
read_in_de <- function(tissue = 'adipose', nominal = F) {
	de = data.frame()
	for (pheno in hist_order) {
		fn = paste0("../../data",
					rep("_liver", (tissue == 'liver')),
					"/topTable_",
					pheno, "_",
					tissue,
					rep("_nominal", nominal),
					".txt")
		newdata = read.table(fn, header = T)
		if (nrow(newdata) > 0) {
			newdata["phenotype"] = pheno
			newdata["gene_ID"] = rownames(newdata)
			rownames(newdata) <- NULL
			de = rbind(de, newdata)
		}
	}
	return(de)
}

# read in significant adipose DE results
de = read_in_de()

# merge on gene name
annot <- read.table("../../data/gencodeV19_annotations_filtered.txt", header = TRUE, sep = '\t')
colnames(annot) <- c("gene_ID", "gene_ID_GENCODEver", "chromosome", "feature", "startPos", "endPos", "strand", "gene_type", "gene_symbol")
annot <- annot[,c("gene_ID", "gene_symbol", "gene_type")]
de = merge(de, annot, by = 'gene_ID')

# remember list of all the DE genes
uniqDEall = unique(de$gene_ID)
write.table(uniqDEall, file = "../../data/allDEgenes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

# remember direction of DE
de['up_in_sick'] = as.character(de$logFC > 0) %>% 
									recode(., 'TRUE' = 'UP', 'FALSE' = 'DOWN')

# add secreted/not secreted column
secreted <- read.table("../../data/protein_class_Secreted_MDSEC.tsv", header = TRUE, sep = "\t", fill = TRUE)
de['secreted'] = de$gene_ID %in% secreted$Ensembl

# add GTEx expression data
gtex <- read.table("../../data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", header = TRUE, skip = 2, sep = '\t')
colnames(gtex)[1] <- "gene_ID_GTExVer"
colnames(gtex)[2] <- "gene_symbol"
novid <- strsplit(gtex$gene_ID_GTExVer, ".", fixed = TRUE)
gtex["gene_ID"] <- sapply(novid, "[[", 1)
# select just subq adipose and liver
gtex = gtex[,c('gene_ID', 'Adipose...Subcutaneous', 'Liver')]
names(gtex) = c('gene_ID', 'adipose_subq', 'liver')
de = merge(de, gtex, by = 'gene_ID', all.x = T)

# remove genes with anything NA at this point
de = de %>% na.omit(.) %>% unique(.)

# add liver DE status
de.liver = read_in_de(tissue = 'liver')
splLiv <- strsplit(de.liver$gene_ID, ".", fixed = TRUE)
de.liver["gene_ID"] <- sapply(splLiv, "[[", 1)
# group on gene ID, aggregate all phenotypes
de.liver = de.liver %>%
						group_by(gene_ID) %>%
						summarise(liver_group = paste(phenotype, collapse = ','))
de['DE_in_liver'] = de$gene_ID %in% de.liver$gene_ID
de = merge(de, de.liver, by = 'gene_ID', all.x = T)

# save table of all DE genes
write.table(de, "../../data/metTable_all_subset.txt", sep = '\t', quote = F, row.names = F)

# select serum biomarker candidates
TPMthresh <- 30
ratioThresh <- 10
de['adipose_expr_pass'] = de$adipose_subq > TPMthresh
de['ratio_thresh_pass'] = de$adipose_subq / de$liver > ratioThresh
# handle divide by 0 error
de[is.na(de$ratio_thresh_pass),'ratio_thresh_pass'] <- T

# identify SBCs based on all 4 filtering steps
de.sbc = de[de$DE_in_liver == F &
						de$secreted == T &
						de$adipose_expr_pass == T &
						de$ratio_thresh_pass == T,]

# remember SBC gene names/IDs
sbc <- unique(de.sbc$gene_symbol)
sbc_ensembl <- unique(de.sbc$gene_ID)
write.table(sbc, file = "../../data/sbc_IDs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(sbc_ensembl, file = "../../data/sbc_IDs_ensembl.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(de.sbc, file = "../../data/metTable_subset_serumBiomarkers.txt", sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

# make a pretty and simple SBC metTable
sbc.clean = de.sbc
sbc.clean['group'] = recode(sbc.clean$phenotype, steatosis = 'Steatosis',
												fibrosis = 'Fibrosis',
												diagnosis = 'NASH')
sbc.clean = sbc.clean %>% 
					group_by(gene_symbol) %>%
					summarise(de_group = paste(group, collapse = ', '),
										secreted = secreted[1],
										DE_in_liver = DE_in_liver[1],
										adipose_median_tpm = adipose_subq[1],
										adipose_over_liver_expr = adipose_subq[1] / liver[1])
colnames(sbc.clean) = c('Gene', 'Adipose DE traits', 'Secreted', 'DE in liver',
												'Adipose expression (median TPM)', 'Adipose:liver expression')
write.table(sbc.clean, file = '../../data/metTable_subset_serumBiomarkers_clean.txt', row.names = F)

# save metTables for all phenotypes
for (p in hist_order) {
	desub = de[de$phenotype == p,]
	write.table(desub, 
			file = paste0("../../data/metTable_", p, "_subset.txt"), 
			sep = '\t', quote = FALSE, row.names = FALSE)
	
	# also save the list of genes for each phenotype if we only care about DE in liver filter
	defilt = desub[desub$DE_in_liver == F,]
	write.table(unique(defilt$gene_ID), 
		paste0("../../data/DE_genes_filter_except_secreted_", p, "_grp.txt"),
		row.names = F, col.names = F, quote = F)
}
# get adipose-aware DE genes across all 3 traits
eefilt_all = de[de$DE_in_liver == F,'gene_ID'] %>% unique(.)
write.table(defilt_all, '../../data/DE_genes_filter_except_secreted_all.txt', row.names = F, col.names = F, quote = F)

# save list of liver-aware DE genes... DE in liver but not adipose
liv.aware = de.liver %>% filter(gene_ID %in% de$gene_ID) %>%
					select(gene_ID) %>% unique()
write.table(liv.aware, '../../data_liver/DE_genes_liver_not_adipose_all.txt', row.names = F, col.names = F, quote = F)




# # plot a venn diagram showing how all the groups fit together
# filtercols = c('DE_in_liver', 'secreted', 
# 					'adipose_expr_pass', 'ratio_thresh_pass')
# de.venn = de[,c('gene_ID', filtercols)] %>% unique(.)
# de.venn['gene_ID'] = NULL

# venncolors = c('gray', 
# 				'dodgerblue1', 
# 				'beige', 
# 				'lightcoral')
# v = ggplot(de.venn, aes(A = DE_in_liver, 
# 						B = secreted, 
# 						C = adipose_expr_pass, 
# 						D = ratio_thresh_pass)) +
# 	# geom_point(aes(x = 0, y = 0), size = 165, color = 'grey') +
# 	geom_venn(show_percentage = F, 
# 						stroke_size = 0.25,
# 						stroke_color = 'white',
# 						fill_color = venncolors,
# 						fill_alpha = 0.8,
# 						set_name_size = 0,
# 						text_size = 5) + 
# 	theme_void() + coord_fixed() + 
# 	coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
# 	theme(plot.margin = margin(c(0, 0, 0, 0), unit = 'in')) +
# 	geom_point(aes(x = 0.5, y = -0.25), size = 25, shape = 1, color = 'black') +
# 	annotate('text', x = 0.5, y = -0.375, label = 'SBCs', fontface = 2) +
# 	annotate('text', x = 0, y = 1, label = as.character(sum(rowSums(de.venn) == 0)), size = 5)

# # add a legend with the colors
# cleannames = c('DE in adipose', 
# 				'DE in liver',
# 				'Secreted',
# 				'Adipose median TPM>30',
# 				'Adipose:liver ratio>10')
# dummy_df = data.frame(foo = factor(cleannames),
# 						bar = 1:5)
# dummy = ggplot(dummy_df, aes(x = bar, y = bar, fill = foo)) +
# 	geom_bar(stat = 'identity') +
# 	scale_fill_manual(breaks = cleannames,
# 						values = c('white', venncolors),
# 						name = '')
# # ggsave('../../figs/foo.png', dummy)
# venn.leg = cowplot::get_legend(dummy)

# # v_out = grid.arrange(v, venn.leg, ncol = 2, nrow = 1, widths = c(10, 1))
# # ggsave("../../figs/sbc_venn.png", v_out, , width = 7.7, height = 7, dpi = 800)

# ggsave("../../figs/sbc_venn.png", v, , width = 7, height = 7, dpi = 800)
# ggsave("../../figs/sbc_venn_legend.png", venn.leg, , width = 2, height = 2, dpi = 800)

# # png('../../figs/sbc_venn2.png', width = 7, height = 7, units = 'in', res = 300)
# # venn(de.venn)
# # dev.off()


