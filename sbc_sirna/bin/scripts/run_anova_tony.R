# run repeated measured ANOVA on Tony's key genes in the control samples

library(edgeR)
library(RNOmni)
library(tidyverse)
library(ggpubr)
library(rstatix)

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

time_order = c('Baseline', '24h', '4D', '7D')

# import cleaned and QCed data
cond = read.table("../../data/condition_per_sample_clean.txt", sep = '\t', header = T)
# gene_info = read.table("../../data/gene_info.txt", sep = '\t', header = T, )
counts = read.table("../../data/knockdown_counts_clean.txt", sep = '\t', header = T)
rownames(counts) = stripVersion(rownames(counts))

# import genes of interest
mygenes = read.table("../../data/tony_genes_mr.txt", header = T)

# select the control samples and interesting genes only
cond = cond[cond$Condition == "Control",]
counts = counts[mygenes$geneID,cond$Sample_ID]

# import the counts and groups into a DGE list
y = DGEList(counts = counts, group = cond$Timepoint)

# normalize by library size
y = calcNormFactors(y)
cpm = cpm(y)

# inverse normal transform the CPMs to meet ANOVA assumptions
cpm_int = data.frame(apply(cpm, 1, RankNorm))
names(cpm_int) = mygenes$gene

# run a repeated measures ANOVA test separately on each gene
anova_results = data.frame()
ttest_results = data.frame()
for (thisgene in mygenes$gene) {
	# gather IDs, timepoints, and response (gene expression) into one table
	expr = cbind(cond[,c("Replicate", "Timepoint")], cpm_int[,thisgene])
	colnames(expr)[3] = thisgene
	expr['Timepoint'] = factor(expr$Timepoint, levels = time_order)

	# compute basic summary stats for each timepoint
	sumstats = expr %>% 
		group_by(Timepoint) %>% 
		get_summary_stats(thisgene, type = 'mean_sd') %>%
		arrange(Timepoint)

	# compute t test across baseline -- 7D
	tt = t.test(expr[expr$Timepoint == "Baseline",thisgene], 
				expr[expr$Timepoint == "7D",thisgene])
	res_tt_new = data.frame(gene = thisgene, pval_baseline_7D = tt$p.value)
	ttest_results = rbind(ttest_results, res_tt_new)

	# plot the data in a boxplot
	boxplot = ggboxplot(expr, x = "Timepoint", y = thisgene, add = 'point') +
		stat_compare_means(comparisons = list(c("Baseline", "7D")), method = 't.test')
	ggsave(paste0('../../figs/', thisgene, '_control_expr_boxplot.png'), boxplot)

	# check for outliers in the data
	expr %>%
		group_by(Timepoint) %>%
		identify_outliers(thisgene)

	# check for normality
	expr %>%
		group_by(Timepoint) %>%
		shapiro_test(thisgene)

	# run the actual ANOVA test
	# res.aov = anova_test(data = expr, dv = thisgene, wid = Replicate, within = Timepoint)
	res.aov = anova_test(data = expr, formula = get(thisgene) ~ Timepoint)
	atable = data.frame(get_anova_table(res.aov))[,c('F', 'p', 'ges')]
	atable['gene'] = thisgene
	anova_results = rbind(anova_results, atable)
}
anova_results['bonf_sig'] = anova_results$p < (0.05 / nrow(mygenes))
anova_results[anova_results$bonf_sig == T,'bonf_sig'] <- '*'
anova_results[anova_results$bonf_sig == F,'bonf_sig'] <- ''

ttest_results['bonf_sig'] = ttest_results$p < (0.05 / nrow(mygenes))
ttest_results[ttest_results$bonf_sig == T,'bonf_sig'] <- '*'
ttest_results[ttest_results$bonf_sig == F,'bonf_sig'] <- ''

write.table(anova_results, '../../data/anova_results_tony_mr.txt', quote = F, row.names = F)
write.table(ttest_results, '../../data/ttest_results_baseline_7D_tony_mr.txt', quote = F, row.names = F)




