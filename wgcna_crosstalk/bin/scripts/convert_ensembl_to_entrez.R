# convert ENSEMBL IDs to entrez IDs
library(biomaRt)

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

# get list of genes I care about
mygenes = unique(c(stripVersion(names(a.datExpr)), stripVersion(names(l.datExpr))))

# query biomart to get the entrez ids
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
filters = listFilters(ensembl)
genes <- getBM(filters="ensembl_gene_id", 
			attributes=c("ensembl_gene_id","entrezgene_id"), 
			values=mygenes, 
			mart=ensembl)

write.table(genes, "../../data/networkgenes_ensembl_to_entrez.txt", quote = F, row.names = F, sep = '\t')
