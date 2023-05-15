library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# prot_order = c('CCDC80')
prot_order = c('CCDC80', 'SOD3')

for (protname in prot_order) {
        # read in hepg2 experiment DE results
        de = read.table(paste0('../../data/topTable_', protname, '.txt'), header = T)

        # clean up the DE table
        de = de %>% 
                # make sure all are significant
                filter(adj.P.Val < 0.05) %>%
                # add a column for protein treatment
                mutate(prot = protname) %>%
                # choose relevant columns
                select(gene_ID, gene_name, logFC, t, 
                                P.Value, adj.P.Val, prot) %>%
                # crop off ensembl ID version
                mutate(gene_ID = stripVersion(gene_ID)) %>%
                # sort on p-value
                arrange(adj.P.Val) %>%
                # convert p-values to scientific notation
                mutate_at(c('P.Value', 'adj.P.Val'), formatC, format = 'e', digits = 3) %>%
                # crop all numbers to 3 decimal places
                mutate(across(where(is.numeric), round, 3)) %>%
                # make the column names look nicer
                rename(`Gene` = gene_name,
                        `Ensembl ID` = gene_ID,
                        `T-statistic` = t,
                        `P-value` = P.Value,
                        `Adj. p-value` = adj.P.Val,
                        `Protein treatment` = prot)
                # output to excel friendly format
                write.csv(de, paste0('../../data/supplement_table_hepg2de_', 
                                                                        protname, '.csv'),
                                                                        row.names = F)

}
