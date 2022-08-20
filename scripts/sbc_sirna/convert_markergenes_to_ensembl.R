library(biomaRt)

# read in adipogenesis marker genes in entrez gene format
adi = read.table("../../data/adipogenesis_genelist.txt", sep = '\t', header = T)

# query biomaRt to get matching ensembl IDs of adipogenesis entrez IDs
# listEnsemblArchives()
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
filters = listFilters(ensembl)
entrezgene = adi$Identifier
genes <- getBM(filters="entrezgene_id", 
			attributes=c("ensembl_gene_id","entrezgene_id"), 
			values=entrezgene, 
			mart=ensembl)

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

write.table(genes$ensembl_gene_id, "../../data/adipogenesis_genelist_ensembl.txt",
			col.names = F, row.names = F, quote = F)

# read in cell type marker genes from Tony
ctm = read.table("../../data/Cryo.snRNA2019.SAT.n15_markers.celltype_padj0.05_uniq.txt", sep = '\t', header = T)

# convert gene symbols to ensembl IDs and remove versions
gene_info = read.table("../../data/gene_info.txt")[,c("gene_name", "gene_type")]
gene_info["gene_id"] = rownames(gene_info)
colnames(gene_info)[1] = "gene"
ctm = merge(ctm, gene_info, by = "gene")

# select preadipocytes and adipocytes
ctm = ctm[ctm$cluster %in% c("Adipocytes", "Preadipocytes"),]

write.table(stripVersion(ctm[ctm$cluster == "Adipocytes",]$gene_id), "../../data/adipocyte_genelist_ensembl.txt",
			col.names = F, row.names = F, quote = F)
write.table(stripVersion(ctm[ctm$cluster == "Preadipocytes",]$gene_id), "../../data/preadipocyte_genelist_ensembl.txt",
                        col.names = F, row.names = F, quote = F)

final_list = c(genes$ensembl_gene_id, stripVersion(ctm$gene_id))

# write out list of ensembl IDs
write.table(final_list, "../../data/adipogenesis_adipocyte_preadipocyte_genelist_ensembl.txt", 
			col.names = F, row.names = F, quote = F)
