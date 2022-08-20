# convert ENSEMBL IDs to entrez IDs
library(biomaRt)

# import lists of genes I care about
sfrp2cor = read.table('../../data/SFRP2_correlated_liver_genes.txt')$V1
allSBCcor = read.table('../../data/allSBC_correlated_liver_genes.txt')$V1
adiDE = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/DE_genes_filter_except_secreted_all.txt')$V1

# download mapping from ensembl to entrez using biomart
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")

labels = c('SFRP2_cor', 'all_SBC_cor', 'adi_aware_DE')
genelists = list(sfrp2cor, allSBCcor, adiDE)

# iterate over each list and convert to entrez ID in MAGENTA-friendly format
for (i in 1:length(labels)) {
	label = labels[i]
	mygenes = genelists[[i]]

	# query biomart object to get the entrez ids
	genes_entrez <- getBM(filters="ensembl_gene_id", 
				attributes=c("ensembl_gene_id","entrezgene_id"), 
				values=mygenes, 
				mart=ensembl)

	# format magenta-style
	magenta = na.omit(genes_entrez$entrezgene_id)
	magenta = c(label, label, magenta)
	magenta = t(magenta)

	write.table(magenta, paste0('../MAGENTA_software_package_vs2_July2011/', label, '_geneset'), 
		row.names = F, col.names = F, quote = F, sep = '\t')
}

