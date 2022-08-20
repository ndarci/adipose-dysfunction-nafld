library(data.table)

# read in expression data for all genes
expr0 = fread("../../data/merged_geneTPM.txt", header = T)

# read in list of sbc genes
sbc = fread("../../data/metTable_subset_serumBiomarkers.txt")$gene_ID

# restrict all expression data to just sbc
expr_sbc = expr0[geneid %in% sbc]

# write result out
fwrite(expr_sbc, "../../data/merged_geneTPM_SBC.txt", sep = '\t')


