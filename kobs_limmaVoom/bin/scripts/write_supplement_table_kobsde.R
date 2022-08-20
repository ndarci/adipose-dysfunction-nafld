library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# import adipose DE results
de.adipose = read.table('../../data/metTable_all_subset.txt', header = T) %>%
		rename(gene_id = gene_ID,
				gene_name = gene_symbol)

# import liver DE results
de.liver = data.frame()
for (pheno in c('steatosis', 'fibrosis', 'diagnosis')) {
	newdata = read.table(paste0('../../data_liver/topTable_', pheno, '_liver.txt'), header = T)
	newdata['phenotype'] = pheno
	de.liver = rbind(de.liver, newdata)
}
de.liver = de.liver %>% mutate(gene_id = stripVersion(rownames(de.liver)))

# read in gene name to ensembl mapping
annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV26_annotations_formatted_genomewide.txt", sep = '\t', header = T)
# add to liver DE information
de.liver = de.liver %>% left_join(annot %>% select(gene_id, gene_name), by = 'gene_id')

# iterate over adipose and liver results
for (tissue in c('adipose', 'liver')) {
	# start with all DE stats
	df = get(paste0('de.', tissue)) %>% 
		# make sure all are significant
		filter(adj.P.Val < 0.05) %>%
		# choose relevant columns
		select(gene_id, gene_name, logFC, t, 
				P.Value, adj.P.Val, phenotype) %>%
		# fix order of phenotypes before sorting
		mutate(phenotype = factor(phenotype, levels = c('steatosis', 'fibrosis', 'diagnosis'))) %>%
		# sort on phenotype, then p-value
		arrange(phenotype, adj.P.Val) %>%
		# make the phenotype values look nicer
		mutate(phenotype = recode(phenotype,
							steatosis = 'Steatosis',
							fibrosis = 'Fibrosis',
							diagnosis = 'NASH')) %>%
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
				`DE phenotype [Steatosis, Fibrosis, or NASH]` = phenotype)
	# output to excel friendly format
	write.csv(df, paste0('../../data/supplement_table_kobsde_', 
								tissue, '.csv'),
								row.names = F)
}






