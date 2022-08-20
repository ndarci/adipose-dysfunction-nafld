library(ggplot2)
library(ggpubr)
library(ggsignif)
library(stringr)
library(tidyverse)
library(edgeR)

theme_set(theme_bw())

# import info about the experiment per sample
cond = read.table("../../data/Illumina_sample_sheet_Uma_Niko_SGBS.csv", sep = ",")
condcols = cond[18,] # get col names
cond = cond[19:nrow(cond),] # remove header
colnames(cond) = condcols

# import covariate data
cov = read.table("../../data/de_cov.txt", header = T)
colnames(cov)[1] = "Sample_ID"
cov$unique_map_pct = as.numeric(cov$unique_map_pct)
cond = merge(cond, cov, by = "Sample_ID") # merge covariates onto condition info

# add individual fields for each experimental condition
cond["Sample_Name"] = str_replace(cond$Sample_Name, "_1_", "_")
condspl = strsplit(cond$Sample_Name, "_", fixed = T)
cond["Timepoint"] = sapply(condspl, "[[", 1)
cond["Condition"] = sapply(condspl, "[[", 2)
cond["Replicate"] = sapply(condspl, "[[", 3)
cond["Timepoint_Condition"] = paste0(cond$Timepoint, "_", cond$Condition)
write.table(cond, "../../data/condition_per_sample_clean.txt", sep = '\t', row.names = F, col.names = T, quote = F)

# import and clean count data from featurecounts
fc_infocols = c("Chr", "Start", "End", "Strand", "Length", "gene_name", "gene_type")
NINFO = length(fc_infocols)
counts = read.table("../../data/featureCounts/gencode19.featureCounts.txt", 
					header = T, row.names = 1, sep = '\t')
splcol = strsplit(colnames(counts)[(NINFO+1):ncol(counts)], ".", fixed = T)
colnames(counts)[(NINFO+1):ncol(counts)] = sapply(splcol, "[[", 4)

# separate out gene metadata
gene_info = counts[,fc_infocols]
counts = counts[,!(colnames(counts) %in% fc_infocols)]
write.table(gene_info, "../../data/gene_info.txt", sep = '\t', row.names = T, col.names = T, quote = F)
write.table(counts, "../../data/knockdown_counts_clean.txt", sep = '\t', row.names = T, col.names = T, quote = F)

# read counts into a DGElist
y = DGEList(counts = counts, group = cond$Timepoint_Condition)

# plot distribution of read counts
readcounts = data.frame(t(y$counts)) %>% 
	select(starts_with("ENSG")) %>% 
	colSums()
rawcountplot = ggplot(mapping = aes(x = readcounts)) + 
	geom_histogram() + 
	scale_y_continuous(trans = "log") +
	xlab("Read counts per gene") + ylab("log(Count)")
ggsave("../../figs/readcount_pergene_distribution_hist.png", rawcountplot)

# get % of total reads each gene represents
readpct = sort(readcounts / sum(readcounts), decreasing = T)
readpct_df = data.frame(gene_ID = names(readpct),
						gene_name = gene_info[names(readpct),]$gene_name,
						prop_of_all_reads = readpct)
head(readpct_df, 10)

# filter out lowly expressed genes
# threshold is smallest number of replicates per condition over total samples when one SBC is excluded

filterdf = data.frame()
for (mincount in seq(0, 20)) {
	keepgenes = filterByExpr(y, min.count = mincount, group = cond$Timepoint_Condition)
	tab = data.frame(table(keepgenes))
	if (TRUE %in% tab$keepgenes) {
		genes_in = tab[tab$keepgenes == T,]$Freq
	} else { genes_in = 0 }
	genes_out = nrow(y) - genes_in
	newdata = data.frame("mincount" = mincount, 
						"num_genes_in" = genes_in,
						"num_genes_out" = genes_out)
	filterdf = rbind(filterdf, newdata)
}
filterplot = ggplot(filterdf, aes(x = mincount, y = num_genes_in)) +
	geom_bar(stat = "identity") +
	xlab(paste0("Minimum count in smallest group size of samples")) +
	ylab(paste0("Number of genes kept (total = ", nrow(y), ")"))
ggsave("../../figs/num_genes_after_filter_bar.png", filterplot)

# just choose default filter for now
keepGenes = filterByExpr(y, group = cond$Timepoint_Condition)
print("Genes kept after expression filter:")
table(keepGenes)
y = y[keepGenes, , keep.lib.sizes = FALSE]

# check gene type of genes we threw out
dropgenes_table = sort(table(gene_info[names(keepGenes[keepGenes == F]),]$gene_type), decreasing = T)
dropgenes_table_pct = round((dropgenes_table / sum(dropgenes_table)), 2)
print("Types of genes dropped by expression filter:")
dropgenes_table_pct

# calculate normalization factors for TMM normalization
y <- calcNormFactors(y)

# compute counts per million to normalize by library size
counts_filt = cpm(y)

# run PCA on the count data
pcaobj = prcomp(t(counts_filt), center = T, scale. = T)

# add first X PCs to the condition data
NPC = 6
firstXpcs = data.frame(pcaobj$x[,1:NPC])
firstXpcs["Sample_ID"] = rownames(firstXpcs)
cond = merge(cond, firstXpcs, by = "Sample_ID")

# check percent variance explained by each PC
pcasum = data.frame(summary(pcaobj)$importance)
print("Variance explained by PCs:")
print(pcasum)

# check covariate correlation with PCs
print("r2 of unique map percentage with PCs:")
rPC = cor(cond$unique_map_pct, firstXpcs[,!(colnames(firstXpcs) == "Sample_ID")], method = "pearson")
r2PC = rPC ^ 2
r2PC

corPCplot = ggplot(mapping = aes(x = colnames(firstXpcs)[1:NPC], y = r2PC)) + 
	geom_bar(stat = "identity") +
	xlab("PC") + ylab("Correlation with Unique Map % (r2)")
ggsave("../../figs/correlation_cov_PCs_bar.png", corPCplot)

# plot PCs and color by experimental condition
time_order = c("Baseline", "24h", "4D", "7D")
cond_order = c("Control", "Controli", "CCDC80", "SOD3")
plotPCs <- function(pcA, pcB) {
	proportions = c()
	for (pc in c(pcA, pcB)) {
		thisprop = round(pcasum["Proportion of Variance",pc]*100, 0)
		proportions = c(proportions, thisprop)
	}

	pcplot = ggplot(cond, aes(x = cond[,pcA], y = cond[,pcB],
					shape = factor(Condition, level = cond_order), 
					colour = factor(Timepoint, level = time_order))) +
			xlab(paste0(pcA, " (", proportions[1], "%)")) +
			ylab(paste0(pcB, " (", proportions[2], "%)")) +
			scale_shape_discrete(name = "Condition") +
			scale_colour_discrete(name = "Timepoint") +
			geom_point(size = 3.5)

	ggsave(paste0("../../figs/", pcA, "_", pcB, ".png"), pcplot, width = 8, height = 7, units = "in")
}

plotPCs("PC1", "PC2")
plotPCs("PC2", "PC3")
plotPCs("PC3", "PC4")
plotPCs("PC4", "PC5")
plotPCs("PC5", "PC6")

# # check if it'll be valid to take mean of replicates (takes a little while to run)
# counts_t = data.frame(t(counts_filt))
# counts_t["Sample_ID"] = rownames(counts_t)
# counts_cond = merge(counts_t, 
# 				cond[,c("Timepoint_Condition", "Replicate", "Sample_ID")], 
# 				by = "Sample_ID")
# counts_sd = counts_cond %>% group_by(Timepoint_Condition) %>%
# 	summarise_at(vars(starts_with("ENSG")), sd)

# # get stats of these SDs
# sdvec = counts_sd %>% select(starts_with("ENSG")) %>% flatten_dbl()
# sdhist_log = ggplot(mapping = aes(x = sdvec)) + 
# 	geom_histogram() +
# 	scale_y_continuous(trans = "log10") +
# 	xlab("SD across replicates") + ylab("log(count)")
# ggsave("../../figs/expression_SD_across_replicates_hist_log.png", sdhist_log)
# sdhist = ggplot(mapping = aes(x = sdvec)) + 
# 	geom_histogram() +
# 	xlab("SD across replicates")
# ggsave("../../figs/expression_SD_across_replicates_hist.png", sdhist)

# check if uniquely mapped % increases with time
ump_barplot = ggbarplot(cond, 
			x = "Timepoint", 
			y = "unique_map_pct",
			add = "mean_sd",
			add_params = list(group = "Timepoint"),
			fill = "lightgrey")
ggsave("../../figs/unique_map_pct_bytimepoint_bar.png", ump_barplot)

# check that the SBC genes were actually knocked down
coolgenes = data.frame("gene_sym" = c("CCDC80", "SOD3", "ADIPOQ", "PPARG"),
						"gene_id" = c("ENSG00000091986.11", "ENSG00000109610.5",
										"ENSG00000181092.5", "ENSG00000132170.15"))
counts_sub = data.frame(t(counts_filt[coolgenes$gene_id,]))
counts_sub["Sample_ID"] = rownames(counts_sub)
cond = merge(cond, counts_sub, by = "Sample_ID")

# measure how much the knockdown genes are expressed relative to control
quantifyKnockdown <- function(cond_sub, gene_sym, gene_id) {
	# remove unimportant columns
	cond_sub = cond_sub[,c("Sample_ID", "Timepoint_Condition", 
						"Timepoint", "Condition", "Replicate", gene_id)]

	# take mean of replicates in same condition
	cond_mean = cond_sub %>% group_by(Timepoint_Condition) %>% 
		summarise(mean_expr = mean(get(gene_id)),
				Timepoint = first(Timepoint),
				Condition = first(Condition))

	# find proportion of knockdown gene expression to control at each time point
	cond_prop = cond_mean %>% 
		pivot_wider(names_from = Condition, values_from = mean_expr) %>%
		group_by(Timepoint) %>%
		summarise(sbc_mean = get(gene_sym)[which(!is.na(get(gene_sym)))],
					# Control_mean = Control[which(!is.na(Control))],
					Controli_mean = Controli[which(!is.na(Controli))]) %>%
		mutate(#sbc_over_control = sbc_mean / Control_mean,
				sbc_over_controli = sbc_mean / Controli_mean,
				gene = gene_sym)

	cond_prop = cond_prop[,c("gene", "Timepoint", 
		# "sbc_over_control",
		"sbc_over_controli")] %>%
		arrange(factor(Timepoint, levels = time_order))

	return(cond_prop)
}

plotSBCknockdown <- function(sbc_sym, gene_id, gene_sym, supplement = FALSE) {
	# remove samples of the other SBC knockdown
	if (supplement == TRUE) {
		# in supplement, include both controls
		cond_sub = cond[cond$Condition %in% c("Control", "Controli", sbc_sym),]
		cond_sub$Condition = factor(cond_sub$Condition, level = c("Control", "Controli", sbc_sym))
		pre = "supplement_"
		legendlabels = c("Control", "Scramble siRNA", paste0(sbc_sym, " siRNA"))
	} else {
		# in main, just keep the scramble control
		cond_sub = cond[cond$Condition %in% c("Controli", sbc_sym),]
		cond_sub$Condition = factor(cond_sub$Condition, level = c("Controli", sbc_sym))
		pre = ""
		legendlabels = c("Scramble siRNA", paste0(sbc_sym, " siRNA"))
	}

	cond_sub$Timepoint = factor(cond_sub$Timepoint, level = time_order)

	# calculate standard deviations of expression for annotation height
	cond_sd = cond_sub %>% 
		group_by(Timepoint_Condition) %>% 
		summarise(sbc_sd = sd(get(gene_id)),
				sbc_mean = mean(get(gene_id)),
				Condition = first(Condition),
				Timepoint = first(Timepoint))
	cond_sd['lowbar'] = cond_sd["sbc_mean"] - cond_sd["sbc_sd"]
	cond_sd['hibar'] = cond_sd["sbc_mean"] + cond_sd["sbc_sd"]
	cond_sd = cond_sd %>% group_by(Timepoint) %>%
		summarise(y.position = max(hibar) + 0.01 * max(hibar))

	# keep track of how tall the plot needs to be
	maxy = max(cond_sd$y.position) + 0.05 * max(cond_sd$y.position)

	# compute t-test statistics
	anno_df = compare_means(as.formula(paste0(gene_id, " ~ Condition")), 
						group.by = "Timepoint",
						method = "t.test",
						ref.group = sbc_sym,
						data = cond_sub[cond_sub$Condition != "Control",]) %>%
			merge(., cond_sd, by = "Timepoint")

	# when looking at SBC expression, add a few extra things to the plot
	if(gene_sym %in% c("CCDC80", "SOD3")) {
		# quantify knockdown efficienty at each time point
		quant_kd = quantifyKnockdown(cond_sub, gene_sym, gene_id)
		# add KD % to the annotation dataframe
		anno_df = anno_df %>%
			merge(., quant_kd, by = "Timepoint") %>%
			mutate(knockdown_pct = round((1 - sbc_over_controli) * 100, 0))
		# include KD % in the annotation label
		plotlabel = "p = {p.format},\n{knockdown_pct}% KD"
	} else {
		plotlabel = "p = {p.format}"
	}

	# generate an expression barplot
	gene_barplot = ggbarplot(cond_sub, 
								x = "Timepoint", 
								y = gene_id,
								fill = "Condition",
								add = "mean_sd", 
								add_params = list(group = "Timepoint_Condition"),
								position = position_dodge(0.8)) +
		xlab("Timepoint of SGBS preadipocyte differentiation") + 
		ylab(paste0(gene_sym, " expression in SGBS cells (CPM)")) +
		coord_cartesian(ylim=c(0, maxy)) +
		scale_fill_discrete(name = "Condition", labels = legendlabels) +
		theme(legend.position = "right") +
		stat_pvalue_manual(anno_df, x = "Timepoint", label = plotlabel, 
				position = position_dodge(0.8))

	ggsave(paste0("../../figs/", pre, sbc_sym, "_knockdown_", gene_sym, "_expression_bar.png"), 
		gene_barplot, width = 9, height = 7, units = "in")
}

plotSBCknockdown("CCDC80", "ENSG00000091986.11", "CCDC80")
plotSBCknockdown("SOD3", "ENSG00000109610.5", "SOD3")
plotSBCknockdown("CCDC80", "ENSG00000091986.11", "CCDC80", T)
plotSBCknockdown("SOD3", "ENSG00000109610.5", "SOD3", T)

# check for evidence of adipogenesis
plotSBCknockdown("CCDC80", "ENSG00000181092.5", "ADIPOQ")
plotSBCknockdown("SOD3", "ENSG00000181092.5", "ADIPOQ")
plotSBCknockdown("CCDC80", "ENSG00000132170.15", "PPARG")
plotSBCknockdown("SOD3", "ENSG00000132170.15", "PPARG")






