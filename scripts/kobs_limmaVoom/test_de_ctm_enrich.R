library(tidyverse)

# import adipose specific DE gene lists
hist_order = c("steatosis_grp", "fibrosis_grp", "diagnosis_grp")#, 
            #"inflammation_grp", "ballooning_grp", "unhealthy_grp")
de = data.frame()
for (hist in hist_order) {
  myGenes <- read.table(paste0("/u/project/pajukant/nikodm/kobs_limmaVoom/data/DE_genes_filter_except_secreted_", hist, ".txt"), header = FALSE)
  myGenes['phenotype'] = hist
  de = rbind(de, myGenes)
}
colnames(de)[1] = 'gene_id'

# import gene name/ID mapping information
annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_formatted_genomewide.txt", sep = '\t', header = T)

# import cell type marker lists
a.marker = read.table('/u/project/pajukant/nikodm/sbc_sirna/data/Cryo.snRNA2019.SAT.n15_markers.celltype_padj0.05_uniq.txt', sep = '\t', header = T)
# add ensembl ID
a.marker['gene_name'] = a.marker$gene
a.marker = merge(a.marker, annot[,c('gene_id', 'gene_name')], by = 'gene_name')

# import list of all genes
allgenes = read.table('../../data/allExprGenes_adipose.txt')$V1
n_allgenes = length(allgenes)

# # do the same for liver module/hepatocytes
# l.marker0 = read.csv('/u/project/pajukant/nikodm/wgcna_crosstalk/data/TableS1.cell_type_markers.csv')
# # select unique cell type marker genes
# l.marker = l.marker0 %>% group_by(Symbol) %>%
#     summarise(Subcell.type = Subcell.type[length(Subcell.type) == 1],
#               Ensembl.ID = Ensembl.ID[length(Subcell.type) == 1],
#               Adjusted.p.value = Adjusted.p.value[length(Subcell.type) == 1]) %>%
#     	rename(cluster = Subcell.type,
#     			gene_name = Symbol,
#     			gene_id = Ensembl.ID,
#     			p_val_adj = Adjusted.p.value)

# test for enrichment
# celltype = 'Adipocytes'
# hist = 'steatosis_grp'
celltypes = unique(a.marker$cluster)

enrich_cols = c('cell_type', 'phenotype', 'n_de', 'n_marker', 'n_overlap', 'pval_enrichment')
enrich_df = data.frame(matrix(NA, nrow = length(celltypes)*length(hist_order), 
						ncol = length(enrich_cols)))
colnames(enrich_df) = enrich_cols

i = 1
for (celltype in celltypes) {
	for (hist in c(hist_order, 'all')) {
		# restrict to current celltype and histology phenotype
		if(hist == 'all') {
			de_sub = de['gene_id'] %>% unique(.)
		} else {
			de_sub = de[de$phenotype == hist,] %>% unique(.)
		}
		a.marker_sub = a.marker[a.marker$cluster == celltype,c('gene_id', 'cluster')] %>% unique(.)

		# count numbers for enrichment test
		nde = nrow(de_sub) # total DE genes
		nctm = nrow(a.marker_sub)# total cell type marker genes
		ndectm = length(intersect(de_sub$gene_id, a.marker_sub$gene_id))# overlap DE-CTM

		# compute enrichment pvalue
		pval_enrich = phyper(ndectm, nctm, n_allgenes-nctm, nde, lower.tail = F)

		enrich_df[i,] = c(celltype, hist, nde, nctm, ndectm, pval_enrich)
		i = i+1
	}
}

enrich_df[enrich_df$pval_enrichment < 0.05,] %>% arrange(pval_enrichment)








