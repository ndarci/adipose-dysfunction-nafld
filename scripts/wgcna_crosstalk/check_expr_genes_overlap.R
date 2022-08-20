stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# load clean count data for each tissue
lnames = load("../../data/adipose_clean_counts_cov.RData")
a.datExpr = datExpr
lnames = load("../../data/liver_clean_counts_cov.RData")
l.datExpr = datExpr

# get the gene lists
a.genes = names(a.datExpr)
l.genes = stripVersion(names(l.datExpr))

# check number of unique and overlapping genes in each list
print(paste('adipose expressed genes:', length(a.genes)))
print(paste('liver expressed genes:', length(l.genes)))
print(paste('genes in common:', length(intersect(a.genes, l.genes))))
print(paste('genes unique to adipose:', sum(!(a.genes %in% l.genes))))
print(paste('genes unique to liver:', sum(!(l.genes %in% a.genes))))
