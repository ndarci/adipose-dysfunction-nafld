library(ggplot2)
library(tidyverse)
library(edgeR)
library(sva)
library(Hmisc)
library(corrplot)

MINCOUNT = 10

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

runDE_timepoint = function(sbc_sym, covcorrect = F, covariate = "", include_all_genes = F) {
	# select samples from this SBC
	cond_sub = cond %>% filter(protein_name == sbc_sym)
	keepsamples = cond_sub$Sample_ID
	y = y0[, keepsamples, keep.lib.sizes = F]

	# filter for expressed genes
	# define gene expression threshold based on number of remaining samples
	keepgenes = filterByExpr(y, min.count = MINCOUNT, group = cond_sub$protein_amt)
	# print(table(keepgenes))
	write.table(stripVersion(rownames(y)[keepgenes]), 
				paste0("../../data/expressed_genes_", 
					sbc_sym, 
					rep("_all_genes_included", include_all_genes), 
					".txt"), 
				sep = '\t', row.names = F, col.names = F, quote = F)
	y = y[keepgenes, , keep.lib.sizes = F]

	cond_sub$protein_amt = factor(cond_sub$protein_amt, level = c('0ng', '20ng'))

	# normalize by library size
	y = calcNormFactors(y)

	# generate the design matrix
	if(covcorrect == TRUE) {
		# # check significance of SV1 calculated for just this timepoint
		# mm = model.matrix(~ Condition, data = cond_sub)
		# nsv = num.sv(cpm(y), mod = mm)
		# print(paste0(tp, "; ", sbc_sym, "; ", 
		# 			rep(paste0("corrected for ", covariate), covcorrect),
		# 			"; significant SVs = ", nsv))
		deformula = as.formula(paste0("~ protein_amt + ", covariate))
	} else {
		deformula = as.formula("~ protein_amt")
	}
	design = model.matrix(deformula, 
			data = cond_sub)
	colnames(design)[2] = "treatment"

	# run the voom-limma pipeline
	v = voom(y, design)
	fit = lmFit(v, design)
	fit2 = eBayes(fit)
	topGenes = topTable(fit2, coef = 2, n = Inf, p.value = Inf, sort.by = "P")

	# check DE gene expression correlation with UM%
	if(dim(topGenes)[1] > 0) {
		de_expr = data.frame(t(cpm(y[rownames(topGenes), , keep.lib.sizes = T])))
		de_expr["Sample_ID"] = rownames(de_expr)
		de_expr = merge(de_expr, cond_sub[,c("Sample_ID", "unique_map_pct")], by = "Sample_ID")
		de_expr_cor = cor(de_expr$unique_map_pct, de_expr[,rownames(topGenes)])
		meandecor = mean(de_expr_cor^2)
	} else {
		meandecor = NA
	}

	return(list(tt = topGenes, meandecor = meandecor))
}

make_volcano_plot <- function(tt, protname) {
	# # create volcano plot given a toptable from limma-voom
	# pthresh <- 0.05
	# tt <- tt %>% mutate(threshold = ifelse(adj.P.Val > pthresh,
	#                                   "Non-significant (adjP>=0.05)",
	#                                   ifelse(logFC > 0, "Up in protein treatment", "Down in protein treatment")))
	# theme_set(theme_bw())
	# volcplot <- ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val))) +
	# 	geom_point(aes(color = threshold), size = 4) +
	# 	scale_color_manual(values = c("Up in protein treatment" = "dodgerblue1", 
	# 		"Down in protein treatment" = "lightcoral", 
	# 		"Non-significant (adjP>=0.05)" = "grey"), 
	# 		breaks = c("Down in protein treatment", "Up in protein treatment", "Non-significant (adjP>=0.05)")) +
	# 	annotate("text", x = 0.81, y = 1.25, label = "Adjusted p-value=0.05", size = 10) +
	# 	geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 1) +
	# 	xlab(paste0("log(FC) in ", protname, " treatment")) + ylab("-log10(Adjusted p-value)") +
	# 	theme(legend.position = c(0.81, 0.9),
	# 		legend.title = element_blank(),
	# 		legend.text = element_text(size = 35),
	# 		text = element_text(family = "Times New Roman"),
	# 		axis.text.x = element_text(size = 30),
	# 		axis.text.y = element_text(size = 30),
	# 		axis.title.x = element_text(size = 40),
	# 		axis.title.y = element_text(size = 40),
	# 		axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
	# 		panel.grid = element_blank())
	# # volcplot
	# fn = paste0("../../figs/volcanoplot_", protname, ".png")
	# ggsave(plot = volcplot, filename = fn, width = 15, height = 13.5, units = "in")
}

# import cleaned data
cond = read.csv('../../data/condition_per_sample_clean.csv') %>%
	rename('Sample_ID' = sample)
cov = read.table("../../data/de_cov.txt", header = T) %>%
	rename('Sample_ID' = sample) %>%
	mutate(unique_map_pct = as.numeric(unique_map_pct))
cond = merge(cond, cov, by = "Sample_ID")
gene_info = read.table("../../data/gene_info.txt", sep = '\t', header = T, )
counts = read.table("../../data/hepg2_counts_clean.txt", sep = '\t', header = T)
targets = read.table('../../data/kobs_liver_DE_hepatocyte_markers_overlap.txt', header = T)$gene_id

# import counts and conditions into a DGElist
y0 = DGEList(counts = counts, group = cond$protein_amt)

# run DE for each timepoint and SBC
# ump_cor_df = data.frame()
for (include_all_genes in c(T, F)) {
	if (include_all_genes == F) {
		# restrict to only overlap of hepatocyte marker genes and liver DE genes from KOBS
		targetgenes_idx = stripVersion(rownames(y0$counts)) %in% targets
		y0 = y0[targetgenes_idx, , keep.lib.sizes = F]
	}
	for (sbc_sym in c("CCDC80", "SOD3")) {
		for (covcorrect in c(T, F)) {
			if (covcorrect == T) {
					covariate = 'unique_map_pct'
				} else {
					covariate = ''
				}
			ext = paste0(sbc_sym, 
						rep(paste0("_covcorrected_", covariate), covcorrect),
						rep("_all_genes_included", include_all_genes))
			# print(ext)
			result = runDE_timepoint(sbc_sym, covcorrect, covariate, include_all_genes)
			# # keep track of correlation of DE gene expression with UM%
			# ump_cor_df = rbind(ump_cor_df, c(sbc_sym, covariate, result$meandecor))
			# add gene symbols
			result = merge(result$tt, gene_info[,c("gene_name", "gene_type")], by = 0)
			colnames(result)[1] = "gene_ID"
			# generate a volcano plot
			make_volcano_plot(result, sbc_sym)
			# save all genes (even non-significant)
			filename = paste0("../../data/topTable_", ext, "_nominal.txt")
			write.table(result, filename, row.names = F, col.names = T, quote = F, sep = '\t')
			# remove non-significant genes after mult. testing correction
			result_filt = result[result$adj.P.Val < 0.05,]
			# save list of DE genes
			delist_filename = paste0("../../data/DE_genes_", ext, ".txt")
			write.table(stripVersion(result_filt$gene_name), delist_filename, 
						row.names = F, col.names = F, quote = F, sep = '\t')
			# save top table
			filename = paste0("../../data/topTable_", ext, ".txt")
			write.table(result_filt, filename, row.names = F, col.names = T, quote = F, sep = '\t')
		}
	}
}

# colnames(ump_cor_df) = c("protein_name", "covariate", "mean_de_ump_cor")
# print(ump_cor_df)

