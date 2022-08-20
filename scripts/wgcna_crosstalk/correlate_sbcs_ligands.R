# library(edgeR)
library(tidyverse)

# add -h to restrict to only healthy people
healthyinput = commandArgs(trailingOnly = T)[1]
healthy = !is.na(healthyinput)
if (healthyinput != '-h') {
	healthy = F 
	print('Unsupported input argument!')
}

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

# # read in adipose expression data
# a.expr0 = read.table("/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt", header = TRUE, row.names = 1)
# a.dge0 = DGEList(a.expr0)

# # read in liver expression data
# l.expr0 = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/gene_counts.match.txt", header = TRUE, row.names = 1)
# l.expr0 <- l.expr0[,7:ncol(l.expr0)]
# l.dge0 = DGEList(l.expr0)

# # restrict to same expressed genes from DE
# a.keepGeneList = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/allExprGenes_adipose.txt")$V1
# a.dge0 = a.dge0[a.keepGeneList, , keep.lib.sizes = F]

# l.keepGeneList = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/allExprGenes_liver.txt")$V1
# l.keepGeneIdx = stripVersion(rownames(l.dge0)) %in% l.keepGeneList
# l.dge0 = l.dge0[l.keepGeneIdx, , keep.lib.sizes = F]

# # calculate CPMs and log-transform them
# a.dge0 = calcNormFactors(a.dge0)
# a.cpm = data.frame(t(cpm(a.dge0, log = T)))

# l.dge0 = calcNormFactors(l.dge0)
# l.cpm = data.frame(t(cpm(l.dge0, log = T)))
# names(l.cpm) = stripVersion(names(l.cpm))

# # get only the individuals in both datasets
# bothindiv = intersect(rownames(a.cpm), rownames(l.cpm))
# a.cpm = a.cpm[bothindiv,]
# l.cpm = l.cpm[bothindiv,]

lnames = load(file = "../../data/adipose_clean_counts_cov.RData")
a.datExpr = datExpr
a.cpm = a.datExpr
a.datTraits = datTraits
lnames = load(file = "../../data/liver_clean_counts_cov.RData")
l.datExpr = datExpr
l.cpm = l.datExpr
names(l.cpm) = stripVersion(names(l.cpm))
l.datTraits = datTraits

if (healthy == T) {
	healthyindiv = a.datTraits %>% 
		filter(steatosis_grp == F & fibrosis_grp == F & diagnosis_grp == F) %>%
		rownames(.)
	a.cpm = a.cpm %>% filter(rownames(a.cpm) %in% healthyindiv)
	l.cpm = l.cpm %>% filter(rownames(l.cpm) %in% healthyindiv)
}

# correlate ADIPOSE expression of all 10 SBCs with 
# LIVER expression of every gene
sbc = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/sbc_IDs_ensembl.txt')$V1

# iterate over SBCs
allcorr_cols = c('r', 'pvalue', 'r2', 'sbc_adipose', 'othergene_liver')
allcorr = matrix(nrow = length(sbc) * ncol(l.cpm), 
		ncol = length(allcorr_cols)) %>%
		data.frame(.)
colnames(allcorr) = allcorr_cols

# function to run each correlation
run_corr <- function(liverexpr, adiexpr) {
	res = cor.test(liverexpr, adiexpr, method = 'pearson')
	return(c(res$estimate, res$p.value))
}

for (i in 1:length(sbc)) {
	id = sbc[i]

	# computate correlations for all genes in liver with this sbc
	corr_res = apply(l.cpm, 2, function(gene) run_corr(gene, a.cpm[,id])) %>% 
		t(.) %>% 
		data.frame(.)
	names(corr_res) = c('r', 'pvalue')

	# add r2
	corr_res['r2'] = (corr_res$r)^2

	# add this gene and correlated gene
	corr_res['sbc_adipose'] = id
	corr_res['othergene_liver'] = rownames(corr_res)

	# concat to big dataframe with other results
	startrow = (i-1) * ncol(l.cpm) + 1
	endrow = startrow + ncol(l.cpm) - 1
	allcorr[startrow:endrow,] = corr_res
}

# add gene symbols
annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_formatted_genomewide.txt", header = TRUE, sep = '\t')
colnames(annot) = c("gene_ID", "gene_ID_GENCODEver", 
	"chromosome", "feature", "startPos", 
	"endPos", "strand", "gene_type", "gene_symbol")
annot = annot[,c("gene_ID", "gene_symbol")]

allcorr = merge(allcorr, annot, by.x = 'sbc_adipose', by.y = 'gene_ID', all.x = T) %>% 
	rename(sbc_adipose_sym = 'gene_symbol') %>%
	merge(., annot, by.x = 'othergene_liver', by.y = 'gene_ID', all.x = T) %>%
	rename(othergene_liver_sym = 'gene_symbol') %>% 
	arrange(r)

write.table(allcorr, paste0('../../data/sbc_all_liver_gene_correlations', rep('_healthy', healthy), '.txt'))












