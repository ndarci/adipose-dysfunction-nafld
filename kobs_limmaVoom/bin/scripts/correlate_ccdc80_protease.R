library(edgeR)
library(tidyverse)
library(Hmisc)

# read in adipose expression data
expr0 = read.table("/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt", header = TRUE, row.names = 1)
dge0 = DGEList(expr0)

# restrict to same expressed genes from DE
keepGeneList = read.table("../../data/allExprGenes_adipose.txt")$V1
dge0 = dge0[keepGeneList, , keep.lib.sizes = F]

# calculate CPMs
dge0 = calcNormFactors(dge0)
cpm = data.frame(t(cpm(dge0, log = F)))

# restrict to protease genes + CCDC80
protease = read.table('../../data/protease_genes.txt')
names(protease) = c('uniprot', 'gene_symbol')
# add CCDC80
protease = rbind(protease, data.frame(uniprot = 'foobar', gene_symbol = 'CCDC80'))

# add gene ID
annot = read.table("../../data/gencodeV19_annotations_formatted_genomewide.txt", header = TRUE, sep = '\t')
colnames(annot) = c("gene_ID", "gene_ID_GENCODEver", "chromosome", "feature", "startPos", "endPos", "strand", "gene_type", "gene_symbol")
annot = annot[,c("gene_ID", "gene_symbol")]

protease = merge(protease, annot, by = 'gene_symbol', all.x = T)

protease = na.omit(protease[protease$gene_ID %in% names(cpm),])

cpm = cpm[,protease$gene_ID]
names(cpm) = protease$gene_symbol

# log2 transform the counts
e <- 1
cpm_log <- data.frame(apply((cpm + e), c(1, 2), log2))

# correlate all these genes with CCDC80
cor <- rcorr(as.matrix(cpm_log), type = "pearson")

cordf = data.frame(r = cor$r[,'CCDC80'], p = cor$P[,'CCDC80']) %>% arrange(r)
cordf = cordf[rownames(cordf) != 'CCDC80',]
cordf['colorcode'] = as.character(as.numeric(cordf$r > 0))
cordf[cordf$p > 0.05,'colorcode'] <- '2'
cordf['gene'] = factor(rownames(cordf), level = rownames(cordf))

print(head(cordf, 10))
print(tail(cordf, 10))

theme_set(theme_bw())
corplot = ggplot(cordf, aes(x = gene, y = r, fill = colorcode)) + 
    scale_fill_manual(name = '',
                        breaks = c('1', '0', '2'),
                        labels = c('+ correlation', '- correlation', 'Non-significant (p>0.05)'),
                        values = c('forestgreen', 'darkorchid2', 'grey')) +
    geom_bar(stat = 'identity') +
    theme(axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1))
ggsave('../../figs/protease_ccdc80_correlations.png', corplot, height = 7, width = 20, dpi = 300)

# check for significant correlations in adipose-aware DE gene list
cordf.sig = cordf %>% filter(p < 0.05/nrow(cordf))

adide = data.frame()
for (pheno in c('steatosis', 'fibrosis', 'diagnosis')) {
	newdata = read.table(paste0('../../data/DE_genes_filter_except_secreted_', pheno, '_grp.txt'))
	newdata['pheno'] = pheno
	colnames(newdata)[1] = 'gene_ID'
	adide = rbind(adide, newdata)
}
adide = merge(adide, annot, by = 'gene_ID', all.x = T)

cordf.sig = merge(cordf.sig, adide, by.x = 'gene', by.y = 'gene_symbol', all.x = T)

print(cordf.sig)




	











