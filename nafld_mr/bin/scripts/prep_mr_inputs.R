library(TwoSampleMR)
library(tidyverse)

# import IV list from LD pruning
iv = read.table('../../data/iv_snps_cond_colocalized_ldpruned.txt')

# overlap with GWAS summary statistics
# tg gwas
gwas.tg = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/TG/TG_bgen.5e-2.txt', header = T)
# nafld gwas
gwas.nafld = read.table('../../data/nafld_gwas_ivonly.txt', header = T)

iv.exp = iv %>% left_join(gwas.tg, by = c('hit1' = 'SNP')) %>%
			select(hit1, BETA, SE, ALLELE1, ALLELE0, A1FREQ,
					CHR, BP, P_BOLT_LMM_INF) %>%
			rename(SNP = 'hit1',
					beta = 'BETA',
					se = 'SE',
					effect_allele = 'ALLELE1',
					other_allele = 'ALLELE0',
					eaf = 'A1FREQ',
					chr = 'CHR',
					position = 'BP',
					pval = 'P_BOLT_LMM_INF') %>%
			mutate(Phenotype = 'TG')

iv.outcome = iv %>% left_join(gwas.nafld, by = c('hit1' = 'SNP')) %>%
			select(hit1, BETA, SE, ALLELE1, ALLELE0, A1FREQ,
					CHR, BP, P_BOLT_LMM_INF) %>%
			rename(SNP = 'hit1',
					beta = 'BETA',
					se = 'SE',
					effect_allele = 'ALLELE1',
					other_allele = 'ALLELE0',
					eaf = 'A1FREQ',
					chr = 'CHR',
					position = 'BP',
					pval = 'P_BOLT_LMM_INF') %>%
			mutate(Phenotype = 'NAFLD')

# convert to TwoSampleMR-friendly format
df.exp = TwoSampleMR::format_data(iv.exp, type = 'exposure')
df.outcome = TwoSampleMR::format_data(iv.outcome, type = 'outcome')

# harmonize effect size directions
df.harm = TwoSampleMR::harmonise_data(exposure_dat = df.exp,
									outcome_dat = df.outcome) %>%
									filter(ambiguous == F)
write.table(df.harm, '../../data/mr_input_table_TwoSampleMR.txt')

# grab some cis and gwas info to highlight coolness of these snps
cis.info = read.table('../../data/iv_snps_cond_colocalized_cis_gwas_stats.txt')
cis.info = cis.info %>% filter(hit1 %in% df.harm$SNP) %>%
                    select(hit1, gene_symbol, beta_eqtl, fdr_eqtl,
                            beta_gwas, p_gwas) %>%
                    unique()
write.table(cis.info, '../../data/final_iv_snps_cis_gwas_stats.txt', row.names = F)






