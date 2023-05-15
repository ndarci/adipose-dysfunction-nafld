library(tidyverse)

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

read_in_kd_de <- function(nominal = F) {
        de = data.frame()
        for (tp in c('Baseline', '24h', '4D', '7D')) {
                for (sbc in c('CCDC80', 'SOD3')) {
                        path = paste0('/u/project/pajukant/nikodm/sbc_sirna/data/topTable_',
                                        tp, '_', 
                                        sbc, 
                                        '_all_genes_included',
                                        rep("_nominal", nominal),
                                        '.txt')
                        newdata = read.table(path, header = T)
                        if (nrow(newdata) > 0) {
                                newdata['timepoint'] = tp
                                newdata['sbc'] = sbc
                                rownames(newdata) <- NULL
                                de = rbind(de, newdata)
                        }
                }
        }
        return(de)
}

# import liver DE data
de.liver = read_in_de(tissue = 'liver')
splLiv <- strsplit(de.liver$gene_ID, ".", fixed = TRUE)
de.liver["gene_ID"] <- sapply(splLiv, "[[", 1)

# merge on gene name
annot <- read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_filtered.txt", header = TRUE, sep = '\t')
colnames(annot) <- c("gene_ID", "gene_ID_GENCODEver", "chromosome", "feature", "startPos", "endPos", "strand", "gene_type", "gene_name")
de.liver = annot %>% select(gene_ID, gene_name) %>%
                right_join(de.liver, by = 'gene_ID')

# import knockdown DE results (with all genes tested)
de.sirna = read_in_kd_de()
# de.sirna = de.sirna[de.sirna$adj.P.Val < 0.05,]

# check for overlaps
ov = de.liver %>%
        # restrict to overlapping genes
        filter(gene_name %in% de.sirna$gene_name) %>% 
        # select sleek subset of columns
        select(gene_name, logFC, adj.P.Val, phenotype) %>%
        # add liver tag to column names
        rename(logFC_liver = logFC,
                adjP_liver= adj.P.Val,
                phenotype_liver = phenotype) %>%
        # merge the sirna data on
        left_join(
                de.sirna %>% 
                select(gene_name, logFC, adj.P.Val, timepoint, sbc) %>%
                rename(logFC_kd = logFC,
                        adjP_kd = adj.P.Val,
                        timepoint_kd = timepoint,
                        sbc_kd = sbc),
                by = 'gene_name'
                ) %>%
        # sort by p value of the knockdown
        arrange(adjP_kd)

# get a more compact version without all the stats
for (sbc in c('CCDC80', 'SOD3')) {
        varname = paste0('ov.', sbc, '.compact')
        assign(varname,
                ov %>% 
                # select this SBC
                filter(sbc_kd == sbc) %>%
                # collapse different DE stats into one line
                group_by(gene_name) %>%
                summarize(DE_phenotypes_liver = paste(unique(phenotype_liver), collapse = ';'),
                        timepoint_kd = paste(unique(timepoint_kd), collapse = ';'),
                        sbc_kd = sbc_kd[1]))
        # add the KD p-values which are the most crucial to check
        assign(varname,
                get(varname) %>%
                left_join(ov %>% select(gene_name, timepoint_kd, sbc_kd, adjP_kd, logFC_kd), by = c('gene_name', 'timepoint_kd', 'sbc_kd')) %>%
                unique() %>%
                mutate(direction_kd = ifelse(logFC_kd > 0, 'up', 'down')) %>%
                arrange(adjP_kd)
        )
        write.table(get(varname), 
                paste0('../../data/potential_liver_target_genes_qPCR_', sbc, '.txt'),
                row.names = F)
}






