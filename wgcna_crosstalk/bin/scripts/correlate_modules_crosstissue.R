library(WGCNA)
library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# import both networks
# adipose
lnames = load(file = "../../data/adipose_clean_counts_cov.RData")
lnames
lnames = load(file = "../../data/adipose_constructed_network.RData")
lnames
a.datExpr = datExpr
a.datTraits = datTraits
a.MEs = MEs 
a.moduleLabels = moduleLabels
a.moduleColors = moduleColors
a.geneTree = geneTree
# liver
lnames = load(file = "../../data/liver_clean_counts_cov.RData")
lnames
lnames = load(file = "../../data/liver_constructed_network.RData")
lnames
l.datExpr = datExpr
l.datTraits = datTraits
l.MEs = MEs 
l.moduleLabels = moduleLabels
l.moduleColors = moduleColors
l.geneTree = geneTree
rm(datExpr, datTraits, MEs, moduleLabels, moduleColors, geneTree)
names(l.datExpr) = stripVersion(names(l.datExpr))

# define x and y dimensions
a.nGenes = ncol(a.datExpr)
l.nGenes = ncol(l.datExpr)
# both tissues have the same number of samples
nSamples = nrow(a.datExpr)

# recalculate MEs with colors
a.MEs = orderMEs(moduleEigengenes(a.datExpr, a.moduleColors)$eigengenes)
l.MEs = orderMEs(moduleEigengenes(l.datExpr, l.moduleColors)$eigengenes)
# keep track of which module comes from which network
colnames(a.MEs) = paste0("A_", colnames(a.MEs))
colnames(l.MEs) = paste0("L_", colnames(l.MEs))

# keep track of how many genes are in each module
get_gene_counts <- function(moduleColors, modNames) {
  me_count = data.frame(table(moduleColors))
  colnames(me_count) = c("module", "count")
  rownames(me_count) = me_count$module
  me_count = me_count[modNames,]
  modNames_count = paste0(modNames, " (", me_count$count, ")")
  return(modNames_count)
}
a.modNames = substring(names(a.MEs), 5)
a.modNames_count = get_gene_counts(a.moduleColors, a.modNames)
l.modNames = substring(names(l.MEs), 5)
l.modNames_count = get_gene_counts(l.moduleColors, l.modNames)

# correlate all MEs across networks
xtissueMEcor = cor(a.MEs, l.MEs, use = "p")
xtissueMEpvalue = corPvalueStudent(xtissueMEcor, nSamples)

# prep for bonferroni correction
a.nModules = dim(xtissueMEpvalue)[1]
l.nModules = dim(xtissueMEpvalue)[2]
adj_p_thresh = 0.05 / (a.nModules * l.nModules)


# # plot heatmap of these correlations
# textMatrix = xtissueMEpvalue < adj_p_thresh
# textMatrix[textMatrix == T] <- ""
# textMatrix[textMatrix == F] <- "X"
# png("../../fig/crosstissue_ME_correlation_heatmap.png", height = 13, width = 13, units = 'in', res = 300)
# par(mar = c(8, 11, 2, 2))
# labeledHeatmap(Matrix = xtissueMEcor,
#               xLabels = paste0("L_ME", l.modNames_count),
#               yLabels = paste0("A_ME", a.modNames_count),
#               colorLabels = FALSE,
#               colors = blueWhiteRed(50),
#               textMatrix = textMatrix,
#               setStdMargins = FALSE,
#               cex.text = 0.5,
#               zlim = c(-1,1))
# dev.off()

# plot custom heatmap with ggplot
library(reshape2)

# calculate first PCs for ordering
a.pc1 = prcomp(xtissueMEcor)$x[,1] %>% 
        data.frame(.) %>%
        rename('pc1' = 1) %>%
        mutate('module' = rownames(.))
l.pc1 = prcomp(t(xtissueMEcor))$x[,1] %>%
        data.frame(.) %>%
        rename('pc1' = 1) %>%
        mutate('module' = rownames(.))

# set up to have gene counts with each module
a.df_mod_count = data.frame(module = rownames(xtissueMEcor),
                            adipose_module_count = a.modNames_count) %>%
                mutate(adipose_module_count = paste0('Adipose ', adipose_module_count)) %>%
                inner_join(a.pc1, by = 'module') %>%
                rename(adipose_module = 'module') %>%
                arrange(pc1)
l.df_mod_count = data.frame(module = names(data.frame(xtissueMEcor)),
                          liver_module_count = l.modNames_count) %>%
                mutate(liver_module_count = paste0('Liver ', liver_module_count)) %>%
                inner_join(l.pc1, by = 'module') %>%
                rename(liver_module = 'module') %>%
                arrange(pc1)

# get data into correct format
            # pivot r values into tidy format
cor.melt = melt(xtissueMEcor) %>%
            # merge with p values
            inner_join(melt(xtissueMEpvalue), by = c('Var1', 'Var2')) %>%
            # clean up column names
            rename(adipose_module = 'Var1',
                    liver_module = 'Var2',
                    r = 'value.x',
                    p = 'value.y') %>%
            # add name plus count column
            inner_join(a.df_mod_count, by = 'adipose_module') %>%
            inner_join(l.df_mod_count, by = 'liver_module') %>%
            # add column for significance annotation
            mutate(nonsig = as.character(p > adj_p_thresh) %>%
                    recode('TRUE' = 'X', 'FALSE' = ''))

# generate the plot
ggcorr = ggplot(cor.melt, aes(x = liver_module_count, y = adipose_module_count, fill = r)) +
    geom_tile() + 
    geom_text(aes(label = nonsig)) +
    # coord_fixed() +
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
    scale_y_discrete(limits = a.df_mod_count$adipose_module_count, position = 'left') +
    scale_x_discrete(limits = l.df_mod_count$liver_module_count) +
    theme(panel.background = element_blank(),
            legend.position = 'right',
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 14))

a.module_order = a.df_mod_count$adipose_module_count
l.module_order = l.df_mod_count$liver_module_count
save(ggcorr, a.module_order, l.module_order, file = '../../fig/gg_crosstissue_ME_correlation_heatmap.RData')
ggsave('../../fig/gg_crosstissue_ME_correlation_heatmap.png', 
    ggcorr, width = 14, height = 14, units = 'in', dpi = 800)







# remember which pairs of modules were significant
get_shared_genes <- function(a.mods, a.all, a.colors, l.mods, l.all, l.colors) {
  out = rep(NA, length(a.mods))
  for (i in seq(1, length(a.mods))) {
    a.genes = a.all[a.colors == a.mods[i]]
    l.genes = l.all[l.colors == l.mods[i]]
    out[i] = length(intersect(a.genes, l.genes))
  }
  return(out)
}
sig.pairs = data.frame(which(xtissueMEpvalue < adj_p_thresh, arr.ind = T))
rownames(sig.pairs) = NULL
sig.pairs = sig.pairs %>% summarise(
                            adipose = rownames(xtissueMEpvalue)[row],
                            liver = colnames(xtissueMEpvalue)[col],
                            r = diag(xtissueMEcor[row, col]),
                            r2 = r^2,
                            p = diag(xtissueMEpvalue[row, col]), 
                            p.adj = p * (a.nModules*l.nModules),
                            adipose_nGenes = table(a.moduleColors)[substring(adipose, 5)],
                            liver_nGenes = table(l.moduleColors)[substring(liver, 5)],
                            sharedGenes = get_shared_genes(substring(adipose, 5), names(a.datExpr), a.moduleColors,
                                                            substring(liver, 5), names(l.datExpr), l.moduleColors)) %>%
                          arrange(p.adj)
write.table(sig.pairs, "../../data/xtissue_ME_correlation_signif.txt", quote = F, row.names = F, sep = '\t')

# compute module membership/connectivity within tissues
# adipose
a.geneModuleMembership = as.data.frame(cor(a.datExpr, a.MEs, use = "p"))
names(a.geneModuleMembership) = paste("MM_A_", a.modNames, sep="")
# compute p values of membership
a.MMPvalue = as.data.frame(corPvalueStudent(as.matrix(a.geneModuleMembership), nSamples))
names(a.MMPvalue) = paste("p.MM_A_", a.modNames, sep="")
# liver
l.geneModuleMembership = as.data.frame(cor(l.datExpr, l.MEs, use = "p"))
names(l.geneModuleMembership) = paste("MM_L_", l.modNames, sep="")
# compute p values of membership
l.MMPvalue = as.data.frame(corPvalueStudent(as.matrix(l.geneModuleMembership), nSamples))
names(l.MMPvalue) = paste("p.MM_L_", l.modNames, sep="")

# save important data for downstream analysis
a.moduleassign = data.frame(gene_id = names(a.datExpr), adipose_module = a.moduleColors) 
l.moduleassign = data.frame(gene_id = names(l.datExpr), liver_module = l.moduleColors) 
save(sig.pairs, a.moduleassign, l.moduleassign, a.geneModuleMembership, l.geneModuleMembership, a.nGenes, l.nGenes, a.modNames, l.modNames, file = '../../data/xtissue_correlation_results.RData') 

