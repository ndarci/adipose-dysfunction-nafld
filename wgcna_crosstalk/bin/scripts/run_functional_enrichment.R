library(anRichment)
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

# recalculate MEs with colors
a.MEs = orderMEs(moduleEigengenes(a.datExpr, a.moduleColors)$eigengenes)
l.MEs = orderMEs(moduleEigengenes(l.datExpr, l.moduleColors)$eigengenes)
# keep track of which module comes from which network
colnames(a.MEs) = paste0("A_", colnames(a.MEs))
colnames(l.MEs) = paste0("L_", colnames(l.MEs))

# import annotations
annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV26_annotations_formatted_genomewide.txt", sep = '\t', header = T)

# import entrez conversion table
conv = read.table("../../data/networkgenes_ensembl_to_entrez.txt", header = T)
a.genes = data.frame(ensembl_gene_id = names(a.datExpr))
a.entrez = merge(a.genes, conv, all.x = T) %>%
			group_by(ensembl_gene_id) %>%
			summarise(entrezgene_id = entrezgene_id[1])
l.genes = data.frame(ensembl_gene_id = stripVersion(names(l.datExpr)))
l.entrez = merge(l.genes, conv, all.x = T) %>%
			group_by(ensembl_gene_id) %>%
			summarise(entrezgene_id = entrezgene_id[1])

# remember cross tissue correlated modules
a.xtissuemodule = "coral1"
l.xtissuemodule = "darkgrey"

# find the GO enrichment of cross-correlated modules
GOcollection = buildGOcollection(organism = "human")
a.GOenr = enrichmentAnalysis(a.moduleColors, a.entrez$entrezgene_id, nBestDataSets = 10, refCollection = GOcollection)
l.GOenr = enrichmentAnalysis(l.moduleColors, l.entrez$entrezgene_id, nBestDataSets = 10, refCollection = GOcollection)


a.GOenr$enrichmentTable[a.GOenr$enrichmentTable$class == a.xtissuemodule,]
l.GOenr$enrichmentTable[l.GOenr$enrichmentTable$class == l.xtissuemodule,]









