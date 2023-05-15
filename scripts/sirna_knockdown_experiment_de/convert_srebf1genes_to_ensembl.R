library(biomaRt)

# read in the pathway gene list
path = read.table('../../data/srebf1_wikipathway.txt', header = T)

# separate out by database
path_ens = path[path$Database == 'Ensembl',]
path_up = path[path$Database == 'Uniprot-TrEMBL',]
path_react = path[path$Database == 'Reactome',]

# convert uniprot IDs to ensembl IDs
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
filters = listFilters(ensembl)
uniprotgene = path_up$Identifier
uniprot_df <- getBM(filters="uniprotswissprot", 
			attributes=c("ensembl_gene_id","uniprotswissprot"), 
			values=uniprotgene, 
			mart=ensembl)

# convert reactome IDs to ensembl IDs
# reactgene = path_react$Identifier
# react_df <- getBM(filters="reactome", 
# 			attributes=c("ensembl_gene_id","reactome"), 
# 			values=reactgene, 
# 			mart=ensembl)

# merge converted data together
outlist = c(path_ens$Identifier, uniprot_df$ensembl_gene_id)

write.table(outlist, "../../data/srebf1_wikipathway_ensembl.txt")