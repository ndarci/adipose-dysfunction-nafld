library(tidyverse)

hist_order = c('steatosis', 'fibrosis', 'diagnosis')

out.table = data.frame()
for (p in hist_order) {
	# import best subsets results
	best = read.table(paste0('../../data/bestSubsets_', p, '.txt'), 
					header = T, sep = '\t')
	best['model'] = p

	# get list of genes considered for this model
	genes = names(best)[1:(which(names(best) == 'r2')-1)]

	# collapse gene membership info into one column
	best['genes_in_model'] = ''
	for (g in genes) {
		best[,g] <- ifelse(best[,g] == '*', paste0(g, '+'), '')
		best['genes_in_model'] = paste0(best[,'genes_in_model'], best[,g])
	}
	best['genes_in_model'] = substr(best$genes_in_model, 1, nchar(best$genes_in_model)-1)

	out.table = best %>% 
			# identify which model was chosen from each phenotype
			mutate(selected = ifelse(best$bic == min(best$bic), '*', '')) %>%
			# choose important columns
			select(model, genes_in_model, r2, adjr2, bic, selected) %>%
			# convert numeric to 3 decimal points
			mutate(across(where(is.numeric), round, 3)) %>%
			# make phenotype labels pretty
			mutate(model = recode(model,
					steatosis = 'Steatosis',
					fibrosis = 'Fibrosis',
					diagnosis = 'NASH')) %>%
			# make the column names pretty
			rename(`Phenotype` = model,
					`Genes in model` = genes_in_model,
					`Adj. r2` = adjr2,
					`BIC` = bic,
					`Best subset` = selected) %>%
			rbind(out.table, .)
}
write.csv(out.table, '../../data/supplement_table_bestsubsets.csv', row.names = F)


