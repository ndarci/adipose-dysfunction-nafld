library(ggplot2)
library(tidyverse)
library(edgeR)
library(sva)
library(Hmisc)
library(corrplot)

MINCOUNT = 10

otherSBC <- function(sbc) {
	if(sbc == "CCDC80") { return("SOD3") }
	else { return("CCDC80") }
}

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

runDE_timepoint = function(tp, sbc_sym, covcorrect = F, covariate = "", include_all_genes = F) {
	# select samples from this timepoint and SBC
	cond_sub = cond[cond$Timepoint == tp,]
	cond_sub = cond_sub[cond_sub$Condition != otherSBC(sbc_sym) &
						cond_sub$Condition != "Control",]
	keepsamples = cond_sub$Sample_ID
	y = y0[, keepsamples, keep.lib.sizes = F]

	# filter for expressed genes
	# define gene expression threshold based on number of remaining samples
	keepgenes = filterByExpr(y, min.count = MINCOUNT, group = cond_sub$Timepoint_Condition)
	# print(table(keepgenes))
	write.table(stripVersion(rownames(y)[keepgenes]), 
				paste0("../../data/expressed_genes_", 
					tp, "_", sbc_sym, 
					rep("_all_genes_included", include_all_genes), 
					".txt"), 
				sep = '\t', row.names = F, col.names = F, quote = F)
	y = y[keepgenes, , keep.lib.sizes = F]

	cond_sub$Condition = factor(cond_sub$Condition, level = c("Controli", sbc_sym))

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
		deformula = as.formula(paste0("~ Condition + ", covariate))
	} else {
		deformula = as.formula("~ Condition")
	}
	design = model.matrix(deformula, 
			data = cond_sub)
	colnames(design)[2] = "knockdown"

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

make_volcano_plot <- function(tt, timepoint, knockdown) {
	# create volcano plot given a toptable from limma-voom
	pthresh <- 0.05
	tt <- tt %>% mutate(threshold = ifelse(adj.P.Val > pthresh,
	                                  "Non-significant (adjP>=0.05)",
	                                  ifelse(logFC > 0, "Up in knockdown", "Down in knockdown")))
	theme_set(theme_bw())
	volcplot <- ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val))) +
		geom_point(aes(color = threshold), size = 4) +
		scale_color_manual(values = c("Up in knockdown" = "dodgerblue1", 
			"Down in knockdown" = "lightcoral", 
			"Non-significant (adjP>=0.05)" = "grey"), 
			breaks = c("Down in knockdown", "Up in knockdown", "Non-significant (adjP>=0.05)")) +
		annotate("text", x = 0.81, y = 1.25, label = "Adjusted p-value=0.05", size = 10) +
		geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 1) +
		xlab(paste0("log(FC) in ", knockdown, " knockdown at ", timepoint)) + ylab("-log10(Adjusted p-value)") +
		theme(legend.position = c(0.81, 0.9),
			legend.title = element_blank(),
			legend.text = element_text(size = 35),
			text = element_text(family = "Times New Roman"),
			axis.text.x = element_text(size = 30),
			axis.text.y = element_text(size = 30),
			axis.title.x = element_text(size = 40),
			axis.title.y = element_text(size = 40),
			axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
			panel.grid = element_blank())
	# volcplot
	fn = paste0("../../figs/volcanoplot_", knockdown, "_", timepoint, ".png")
	ggsave(plot = volcplot, filename = fn, width = 15, height = 13.5, units = "in")
}

# import cleaned data
cond = read.table("../../data/condition_per_sample_clean.txt", sep = '\t', header = T)
gene_info = read.table("../../data/gene_info.txt", sep = '\t', header = T, )
counts = read.table("../../data/knockdown_counts_clean.txt", sep = '\t', header = T)
adi = read.table("../../data/adipogenesis_adipocyte_preadipocyte_genelist_ensembl.txt", header = F)$V1

# import counts and conditions into a DGElist
y0 = DGEList(counts = counts, group = cond$Condition)

# # compute SVs of the full count matrix (separated by SBC)
# for (sbc in c("CCDC80", "SOD3")) {
# 	# exclude samples from other SBC and regular control
# 	cond_sub = cond[cond$Condition != otherSBC(sbc) &
# 						cond$Condition != "Control",]
# 	keepsamples = cond_sub$Sample_ID
# 	y_sva = y0[, keepsamples, keep.lib.sizes = F]
# 	cond_sub$Condition = factor(cond_sub$Condition, level = c("Controli", sbc))

# 	# exclude lowly expressed genes
# 	keepgenes = filterByExpr(y_sva, min.count = MINCOUNT, group = cond_sub$Timepoint_Condition)
# 	y_sva = y_sva[keepgenes, , keep.lib.sizes = F]

# 	# calculate normalization factors for CPM
# 	y_sva <- calcNormFactors(y_sva)

# 	# compute SVs
# 	mm = model.matrix(~ Condition + Timepoint, data = cond_sub)
# 	svaobj = sva(cpm(y_sva), mod = mm)

# 	# compute PCs
# 	pcaobj = prcomp(t(cpm(y_sva)), center = T, scale. = T)

# 	# wrap SV1 and sample IDs into a dataframe
# 	sv_df = data.frame(cbind(cond_sub$Sample_ID, svaobj$sv[,1], pcaobj$x[,1]))
# 	colnames(sv_df) = c("Sample_ID", paste0("SV1_", sbc), paste0("PC1_", sbc))

# 	# merge library sizes, SV1, and PC1 onto the sample info dataframe
# 	libsize = y_sva$samples 
# 	libsize["Sample_ID"] = rownames(libsize) 
# 	cond_sub = merge(cond_sub, sv_df, by = "Sample_ID") %>% 
# 		merge(., libsize[,c("Sample_ID", "lib.size")], by = "Sample_ID")

# 	# correlate SV1 with desired tech factors
# 	cor_df = cond_sub[,c(paste0("SV1_", sbc), paste0("PC1_", sbc), 
# 						"unique_map_pct", "lib.size")]
# 	cor = rcorr(as.matrix(cor_df), type = "pearson")
# 	bonfsig = 0.05/(ncol(cor$r))
# 	cor$P[which(is.na(cor$P))] <- 0
# 	png(paste0("../../figs/tech_factor_corrplot_", sbc, ".png"), width = 2300, height = 2300, res = 300)
# 	corrplot.mixed(cor$r,
# 		p.mat = cor$P, sig.level = bonfsig,
# 		pch.cex = 5, pch.col = "gray",
# 		tl.cex = 1, tl.col = "black",
# 		cl.cex = 1,
# 		lower = "number",
# 		upper = "color")
# 	dev.off()

# 	# merge this new df onto cond (Control and other SBC will be NA)
# 	cond = merge(cond, sv_df, by = "Sample_ID", all.x = T)
# 	cond[,paste0("SV1_", sbc)] = as.numeric(cond[,paste0("SV1_", sbc)])
# 	cond[,paste0("PC1_", sbc)] = as.numeric(cond[,paste0("PC1_", sbc)])
# }

# run DE for each timepoint and SBC
ump_cor_df = data.frame()
for (include_all_genes in c(T, F)) {
	if (include_all_genes == F) {
		# restrict to only the adipogenesis marker genes (and knockdown genes)
		adigenes_idx = stripVersion(rownames(y0$counts)) %in% c(adi, "ENSG00000091986", "ENSG00000109610")
		y0 = y0[adigenes_idx, , keep.lib.sizes = F]
	}
	for (tp in c("Baseline", "24h", "4D", "7D")) {
		for (sbc_sym in c("CCDC80", "SOD3")) {
			# for (covariate in c(paste0("SV1_", sbc_sym), "unique_map_pct", "")) {
			# 	if (covariate == "") {
			# 		covcorrect = F 
			# 	} else { covcorrect = T }
				covariate = ""
				covcorrect = F
				ext = paste0(tp, "_", 
							sbc_sym, 
							rep(paste0("_covcorrected", covariate), covcorrect),
							rep("_all_genes_included", include_all_genes))
				# print(ext)
				result = runDE_timepoint(tp, sbc_sym, covcorrect, covariate, include_all_genes)
				# keep track of correlation of DE gene expression with UM%
				ump_cor_df = rbind(ump_cor_df, c(tp, sbc_sym, covariate, result$meandecor))
				# add gene symbols
				result = merge(result$tt, gene_info[,c("gene_name", "gene_type")], by = 0)
				colnames(result)[1] = "gene_ID"
				# generate a volcano plot
				make_volcano_plot(result, tp, sbc_sym)
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
			# }
		}
	}
}

colnames(ump_cor_df) = c("timepoint", "knockdown", "covariate", "mean_de_ump_cor")
# print(ump_cor_df)

