sex = commandArgs(trailingOnly = T)[1]

# read in de results
de = read.table(paste0("../../data_liver/topTable_steatosis_", sex, "_liver.txt"), header = T)

# add up in cases column
de['up_in_cases'] = de$logFC > 0
de[de$up_in_cases == T,'up_in_cases'] <- "UP"
de[de$up_in_cases == F,'up_in_cases'] <- "DOWN"

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# add gene symbol
annot = read.table("../../data/gencodeV26_annotations_formatted_genomewide.txt", header = T)
de['gene_id'] = stripVersion(rownames(de))
de = merge(de, annot[,c('gene_id', 'gene_name', 'gene_type')], by = 'gene_id', all.x = T)

# write a nice table out
de.out = de[,c('gene_name', 'gene_id', 'up_in_cases', 'logFC', 'adj.P.Val', 'gene_type')]
write.table(de.out, paste0("../../data_liver/topTable_steatosis_", sex, "_liver_formatted.txt"), sep = '\t', quote = F, row.names = F)
