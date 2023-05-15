library(ggplot2)
library(ggpubr)
library(ggsignif)
library(stringr)
library(tidyverse)
library(edgeR)

theme_set(theme_bw())

amt_order = c('0ng', '20ng')
prot_order = c("CCDC80", "SOD3")

# import info about the experiment per sample
cond = read.csv('../../data/condition_per_sample_clean.csv') %>%
	rename('Sample_ID' = sample) %>%
	mutate(prot_name_amt = paste(protein_name, protein_amt, sep = '_'))

# import covariate data
cov = read.table("../../data/de_cov.txt", header = T) %>%
	rename('Sample_ID' = sample) %>%
	mutate(unique_map_pct = as.numeric(unique_map_pct))
cond = merge(cond, cov, by = "Sample_ID") # merge covariates onto condition info

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
write.table(counts, "../../data/hepg2_counts_clean.txt", sep = '\t', row.names = T, col.names = T, quote = F)

# read counts into a DGElist
y = DGEList(counts = counts, group = cond$prot_name_amt)

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
	keepgenes = filterByExpr(y, min.count = mincount, group = cond$prot_name_amt)
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
keepGenes = filterByExpr(y, group = cond$prot_name_amt)
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

# check how much the unique map % changes across experimental condition
# do separately for CCDC80 and SOD3
for (protname in prot_order) {
	cond.sub = cond %>% filter(protein_name == protname)
	ump_stat = cond.sub %>%
			group_by(protein_amt) %>%
			summarise(mean_unique_map_pct = mean(unique_map_pct),
						sd_unique_map_pct = sd(unique_map_pct))
	print(paste0("T-test between unique map percent in treatment and control groups, ", protname))
	print(ump_stat)
	print(t.test(cond.sub[cond.sub$protein_amt == '0ng',]$unique_map_pct, 
			cond.sub[cond.sub$protein_amt == '20ng',]$unique_map_pct))
}

# plot PCs and color by experimental condition
plotPCs <- function(pcA, pcB) {
	proportions = c()
	for (pc in c(pcA, pcB)) {
		thisprop = round(pcasum["Proportion of Variance",pc]*100, 0)
		proportions = c(proportions, thisprop)
	}

	pcplot = ggplot(cond, aes(x = cond[,pcA], y = cond[,pcB],
					shape = factor(protein_amt, level = amt_order),
					color = factor(protein_name, level = prot_order))) +
			xlab(paste0(pcA, " (", proportions[1], "%)")) +
			ylab(paste0(pcB, " (", proportions[2], "%)")) +
			scale_shape_discrete(name = "Protein amount") +
			scale_color_discrete(name = "Protein name") +
			geom_point(size = 3.5)

	ggsave(paste0("../../figs/", pcA, "_", pcB, ".png"), pcplot, width = 8, height = 7, units = "in")
}

plotPCs("PC1", "PC2")
plotPCs("PC2", "PC3")
plotPCs("PC3", "PC4")
plotPCs("PC4", "PC5")
plotPCs("PC5", "PC6")






