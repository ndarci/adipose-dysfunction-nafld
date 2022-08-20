library(tidyverse)

time_order = c('Baseline', '24h', '4D', '7D')

# read in DE data per timepoint/sbc
# do two tables, one per knockdown
read_in_de = function(kd) {
	df = data.frame()
	for (timepoint in time_order) {
		newdata = read.table(paste0('../../data/metTable_', timepoint, '_', kd, '.txt'), header = T, sep = '\t')
		df = rbind(df, newdata)
	}
	return(df)
}

# iterate over each sbc
for (sbc in c('CCDC80', 'SOD3')) {
	# get all the DE data from this knockdown
	de = read_in_de(sbc) %>%
		rename(gene_id = gene_ID) %>%
		# make sure all are significant
		filter(adj.P.Val < 0.05) %>%
		# choose relevant columns
		select(gene_id, gene_name, logFC, t, 
				P.Value, adj.P.Val, timepoint) %>%
		# fix order of phenotypes before sorting
		mutate(timepoint = factor(timepoint, levels = time_order)) %>%
		# sort on phenotype, then p-value
		arrange(timepoint, adj.P.Val) %>%
		# convert p-values to scientific notation
		mutate_at(c('P.Value', 'adj.P.Val'), formatC, format = 'e', digits = 3) %>%
		# crop all numbers to 3 decimal places
		mutate(across(where(is.numeric), round, 3)) %>%
		# make the column names look nicer
		rename(`Gene` = gene_name,
				`Ensembl ID` = gene_id,
				`T-statistic` = t,
				`P-value` = P.Value,
				`Adj. p-value` = adj.P.Val,
				`DE timepoint [Baseline, 24h, 4D, or 7D]` = timepoint)
	# write out in an excel friendly format
	write.csv(de, paste0('../../data/supplement_table_knockdownde_', sbc, '.csv'),
					row.names = F)
}










