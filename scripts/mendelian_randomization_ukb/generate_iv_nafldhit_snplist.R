library(tidyverse)

# import candidate IVs
iv = read.table('../../data/cond_coloc_IVregion_cis-eQTL_UKBgwasTG_0.8PPH4.txt')

# get genetic positions of the IV candidates
iv.pos = data.frame()
for (gene in iv$gene_ID %>% unique()) {
	newdata = read.table(paste0('../../data/snp_overlap/TG_noNAFLD_', gene, '_SNPoverlap.txt'))
	iv.pos = rbind(iv.pos, newdata)
}
iv.pos = iv.pos %>% filter(rsID_gwas %in% iv$hit1) %>%
				# select(rsID_gwas, a1_gwas, a0_gwas, chr_eqtl, pos_eqtl) %>%
				unique()

iv = iv %>% left_join(iv.pos, by = c(c('hit1' = 'rsID_gwas'), c('gene_ID' = 'gene_eqtl')))
write.table(iv, '../../data/iv_snps_cond_colocalized_cis_gwas_stats.txt')

# import nafld gwas results
gwas.nafld = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen_5e-2.txt', header = T)

# select significant gwas signals
gwas.nafld = gwas.nafld %>% filter(P_BOLT_LMM_INF < 5e-8)

# select variables needed for snplist (METSIM SNP ID format)
iv = iv %>% select(hit1, a1_gwas, a0_gwas) %>% unique()
gwas.nafld = gwas.nafld %>% select(SNP, ALLELE1, ALLELE0)
names(gwas.nafld) = names(iv)

# create the table which will be used to generate the output list 
snplist = rbind(iv, gwas.nafld)

write.table(iv, '../../data/iv_snps_cond_colocalized_clean.txt')

write.table(snplist, '../../data/iv_candidate_snps_plus_nafld_gwashits.txt', sep = ':',
			col.names = F, row.names = F, quote = F)

