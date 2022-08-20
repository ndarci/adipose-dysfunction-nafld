## runSBCcorrelation.R
# run a gene x gene association test for all the SBCs and their expression values

library(ggplot2)
library(Hmisc)
# library(corrplot)
library(edgeR)

# read in adipose expression data
expr0 = read.table("/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt", header = TRUE, row.names = 1)
dge0 = DGEList(expr0)

# restrict to same expressed genes from DE
keepGeneList = read.table("../../data/allExprGenes_adipose.txt")$V1
dge0 = dge0[keepGeneList, , keep.lib.sizes = F]

# calculate CPMs
dge0 = calcNormFactors(dge0)
tmm = data.frame(t(cpm(dge0, log = F)))

# restrict to SBCs
sbc = data.frame(gene_id = read.table("../../data/sbc_IDs_ensembl.txt")$V1,
                gene_name = read.table("../../data/sbc_IDs.txt")$V1)
tmm = tmm[,sbc$gene_id]
names(tmm) = sbc$gene_name

# log2 transform the counts
e <- 1
tmm_log <- data.frame(apply((tmm + e), c(1, 2), log2))

# # plot the counts of 5 random genes to see if they're normally distributed
# ggplot(data = tmm, aes(x = CD300LG)) + geom_histogram()
# ggplot(data = tmm, aes(x = VEGFB)) + geom_histogram()
# ggplot(data = tmm, aes(x = TIMP3)) + geom_histogram()
# #ggplot(data = tmm, aes(x = DPT)) + geom_histogram()
# #ggplot(data = tmm, aes(x = MGP)) + geom_histogram()
# # ... they're not
# # same with log-counts
# ggplot(data = tmm_log, aes(x = CD300LG)) + geom_histogram()
# ggplot(data = tmm_log, aes(x = VEGFB)) + geom_histogram()
# ggplot(data = tmm_log, aes(x = TIMP3)) + geom_histogram()
# #ggplot(data = tmm_log, aes(x = DPT)) + geom_histogram()
# #ggplot(data = tmm_log, aes(x = MGP)) + geom_histogram()
# # ... those look better

# compute the correlations
cor <- rcorr(as.matrix(tmm_log), type = "pearson")
# save each component
write.table(cor$r, "../../data/geneByGeneCorrSBC_r.txt", quote = F)
write.table(cor$P, "../../data/geneByGeneCorrSBC_p.txt", quote = F)
write.table(cor$n, "../../data/geneByGeneCorrSBC_n.txt", quote = F)

# find independently expressed genes
corthresh <- 0.5
isLoCor <- data.frame(apply(cor$r, c(1, 2), function (x) x < corthresh && x > -1*corthresh))
isLoCor["independent"] <- rowSums(isLoCor) == (ncol(isLoCor) - 1)
indepgenes <- rownames(isLoCor[isLoCor$independent == T,])
indepgenes

# find the most co-expressed gene
# add the absolute correlation coefficients of all the genes 
# (include self-correlation b/c it affects them all the same)
sumcor <- rowSums(abs(cor$r))
sumcor <- sort(sumcor, decreasing = T)
sumcor2 <- rowSums((cor$r)^2)
sumcor2 <- sort(sumcor2, decreasing = T)
# what's the most co-expressed gene?
sumcor[1]
sumcor2[1]

# fill in NA values
for (i in seq(1, nrow(cor$P))) {
cor$P[i,i] = 0
}

# # plot correlations
# scaledodgercoral = colorRampPalette(c('lightcoral', 'white', 'dodgerblue1'))(10)
# png("../../figs/geneByGeneCorrSBC.png", width = 7, height = 7, units = 'in', res = 800)
# corrplot(cor$r, type = "full", p.mat = cor$P, sig.level = 0.05/nrow(cor$r),
#          method = "color", order = "FPC", 
#          col = scaledodgercoral,
#          pch.cex = 1.75, pch.col = "black",
#          tl.cex = 1, tl.col = "black")
# # colorlegend(scaledodgercoral, c(seq(-1, 1, 0.2)), 
# #             align = 'c', vertical = F, addlabels = T)
# dev.off()

# try it custom with ggplot
library(reshape2)
library(tidyverse)
library(cowplot)
# library(gridExtra)

# calculate first pc for ordering
pc1 = prcomp(cor$r)$x[,1] %>% sort(.)
pc1_ord = names(pc1)

# get data into correct format
cor.melt = merge(melt(cor$r), melt(cor$P), by = c('Var1', 'Var2'))
names(cor.melt) = c('gene1', 'gene2', 'r', 'p')
cor.melt['nonsig'] = as.character(cor.melt$p > 0.05/nrow(cor$r)) %>% 
                        recode(., 'TRUE' = 'X',
                                'FALSE' = '')

# generate the plot
ggcorr = ggplot(cor.melt, aes(x = gene1, y = gene2, fill = r)) +
    geom_tile() + 
    geom_text(aes(label = nonsig)) +
    coord_fixed() +
    scale_fill_gradient2(low = 'darkorchid2', 
                        mid = 'white', 
                        high = 'forestgreen',
                        limits = c(-1, 1)) +
    guides(fill = guide_colorbar(title = expression(italic('R')),
                                barwidth = 0.5, 
                                barheight = 15,
                                label.position = 'right',
                                title.hjust = 0,
                                ticks = F)) + 
    scale_y_discrete(limits = pc1_ord, position = 'left') +
    scale_x_discrete(limits = pc1_ord) +
    theme(panel.background = element_blank(),
            legend.position = 'right',
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 14))

ggsave('../../figs/ggcorr_sbc.png', ggcorr, width = 7.7, height = 7, units = 'in', dpi = 800)

# # get legend
# ggcor.leg = cowplot::get_legend(ggcorr)

# out = grid.arrange(ggcor.leg, ggcorr + theme(legend.position = 'none'),
#                     nrow = 1, ncol = 2, widths = c(1, 10))

# ggsave('../../figs/ggcorr_sbc.png', out, width = 7.7, height = 7, units = 'in', dpi = 200)















