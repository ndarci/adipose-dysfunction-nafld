# read in all the magenta outputs for each gwas/gene set combo
mag = data.frame()
for (gwas in c('TG', 'NAFLD_imp', 'ALT', 'GGT')) {
	for (geneset in c('SFRP2_cor', 'all_SBC_cor', 'adi_aware_DE')) {
		resfile = paste0('../MAGENTA_software_package_vs2_July2011/',
						'Output_MAGENTA_',
						gwas, '_gwas_',
						geneset, '_geneset_10000perm_Jul25_22/',
						'MAGENTA_pval_GeneSetEnrichAnalysis_',
						gwas, '_gwas_',
						geneset, '_geneset_1000kb_upstr_1000kb_downstr_10000perm_Jul25_22.results')
		newdata = read.table(resfile, header = T, comment.char = '', fill = NA)
		newdata['gwas'] = gwas
		newdata['geneset'] = geneset
		mag = rbind(mag, newdata)
	}
}

write.table(mag, '../../data/magenta_results_all.txt', row.names = F)