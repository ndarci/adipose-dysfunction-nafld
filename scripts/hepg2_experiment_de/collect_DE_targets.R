library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# function to read in kobs DE results
hist_order = c('steatosis', 'fibrosis', 'diagnosis')
read_in_de <- function(tissue = 'adipose', nominal = F) {
        de = data.frame()
        for (pheno in hist_order) {
                fn = paste0("/u/project/pajukant/nikodm/kobs_limmaVoom/data",
                                        rep("_liver", (tissue == 'liver')),
                                        "/topTable_",
                                        pheno, "_",
                                        tissue,
                                        rep("_nominal", nominal),
                                        ".txt")
                newdata = read.table(fn, header = T)
                if (nrow(newdata) > 0) {
                        newdata["phenotype"] = pheno
                        newdata["gene_ID"] = rownames(newdata)
                        rownames(newdata) <- NULL
                        de = rbind(de, newdata)
                }
        }
        return(de)
}

# read in the liver DE results from KOBS
de.liver = read_in_de(tissue = 'liver') %>%
                mutate(gene_id = stripVersion(gene_ID))

# read in the adipose DE results from KOBS, just to see where CCDC80/SOD3 were DE
de.adipose = read_in_de() 
name_table = data.frame(gene_ID = c('ENSG00000109610', 'ENSG00000091986'), 
                        gene_name = c('SOD3', 'CCDC80'))
de.adipose = de.adipose %>%
                left_join(name_table, by = 'gene_ID') %>%
                filter(gene_name %in% c('SOD3', 'CCDC80'))
# OK, they're both DE in all 3 of them

# read in the liver cell type marker gene data (unique CTM genes only)
l.marker = read.table('/u/project/pajukant/nikodm/wgcna_crosstalk/data/seur.CellType.markers.unique.txt', header = T)
l.marker = l.marker %>% 
                mutate(gene_id = stripVersion(gene)) %>%
                # select only hepatocyte markers
                filter(cluster %in% paste0('Hep-', 1:14))

# overlap liver DE genes with hepatocyte CTM genes
de.ctm.ov = de.liver %>%
                select(gene_id, logFC, P.Value, adj.P.Val, phenotype) %>%
                rename(logFC_de = logFC,
                        pval_de = P.Value,
                        pval_adj_de = adj.P.Val,
                        phenotype_de = phenotype) %>%
                inner_join(l.marker %>% 
                                select(gene_id, gene_name, 
                                        cluster, p_val, p_val_adj) %>%
                                rename(cluster_ctm = cluster,
                                        pval_ctm = p_val,
                                        pval_adj_ctm = p_val_adj), 
                                by = 'gene_id')

# count the target spaces we'd be talking about
n_liver_de = length(unique(de.liver$gene_id))
n_hep_ctm = length(unique(l.marker$gene_id))
n_liver_de_and_hep_ctm = length(unique(de.ctm.ov$gene_id))

print(paste('Genes DE in liver in KOBS (all NAFLD phenotypes):', n_liver_de))
print(paste('Hepatocyte marker genes (all subtypes):', n_hep_ctm))
print(paste('Genes shared between these two lists:', n_liver_de_and_hep_ctm))

write.table(de.ctm.ov, '../../data/kobs_liver_DE_hepatocyte_markers_overlap.txt', row.names = F)



