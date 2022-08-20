library(tidyverse)

# import sbc-liver gene correlation results
c = read.table('../../data/xtissue_ligand_cor_stats.txt', header = T)

# select only the robust correlations
ch = c %>% filter(nom_sig_in_healthy == T)

# check that the direction of correlation is preserved in healthy indiv
print('direction of correlation preserved:')
all((ch$r > 0) == (ch$healthy_r > 0))

# count number of correlations per sbc
ch.sbc = ch %>% group_by(sbc_adipose_sym) %>%
				summarise(ncor = length(othergene_liver)) %>%
				arrange(ncor)
print('liver correlations per sbc:')
ch.sbc

# count number of liver DE genes
ch.liverde = ch %>% select(othergene_liver, starts_with('liverDE_adjP')) %>%
					unique %>%
					apply(., 2, function(col) !is.na(col)) %>% colSums
print('liver DE genes per phenotype:')
ch.liverde

# check for interesting wgcna modules
cooladipose = c('lightyellow', 'cyan')
coolliver = c('saddlebrown')
ch.wgcna = ch %>% filter(sbc_adipose_module %in% cooladipose | 
					othergene_liver_module %in% coolliver) %>%
				select(sbc_adipose_sym, othergene_liver_sym, r, pvalue, r2, 
						sbc_adipose_module, othergene_liver_module)
print('correlated genes in cool wgcna modules:')
ch.wgcna

# count cis gwas signals 
ch.gwas = ch %>% select(starts_with('cis_gwas')) %>% 
				na.omit %>% 
				colSums
print('gwas hits per phenotype:')
ch.gwas

# check for ccdc80 and sod3 correlations
ch.kd = ch %>% filter(sbc_adipose_sym %in% c('CCDC80', 'SOD3'))
print('CCDC80 and SOD3 correlations:')
ch.kd

# output gene lists for pathway enrichment
write.table(unique(ch$othergene_liver), '../../data/allSBC_correlated_liver_genes_healthy.txt', row.names = F, col.names = F, quote = F)
write.table(unique(ch %>% filter(sbc_adipose_sym == 'SFRP2') %>% .$othergene_liver), '../../data/SFRP2_correlated_liver_genes_healthy.txt', row.names = F, col.names = F, quote = F)




