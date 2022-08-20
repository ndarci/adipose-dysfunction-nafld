library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# liver
l.marker0 = read.csv('/u/project/pajukant/nikodm/wgcna_crosstalk/data/TableS1.cell_type_markers.csv')
# select unique cell type marker genes
l.marker = l.marker0 %>% group_by(Symbol) %>%
    summarise(Subcell.type = Subcell.type[length(Subcell.type) == 1],
              Ensembl.ID = Ensembl.ID[length(Subcell.type) == 1],
              Adjusted.p.value = Adjusted.p.value[length(Subcell.type) == 1]) %>%
      rename(cluster = Subcell.type,
          gene_name = Symbol,
          gene_id = Ensembl.ID,
          p_val_adj = Adjusted.p.value)
l.marker['gene_id'] = stripVersion(l.marker$gene_id)

for (celltype in paste0('Hep_', 1:3)) {
  genelist = l.marker[l.marker$cluster == celltype,]$gene_id
  write.table(genelist, paste0('../../data/liver_ctm_genes_', celltype, '.txt'), quote = F, row.names = F, col.names = F)
}