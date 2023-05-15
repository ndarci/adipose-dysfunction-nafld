library(tidyverse)

# read in adipose and liver network stats
lnames = load('../../data/xtissue_correlation_results.RData')

# read in gene name to ensembl mapping
annot = read.table("/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV26_annotations_formatted_genomewide.txt", sep = '\t', header = T)

# iterate over all cool modules
for (thing in list(c('Adipose', 'lightyellow'), c('Adipose', 'cyan'),
				c('Liver', 'saddlebrown'), c('Adipose', 'lavenderblush3'),
				c('Adipose', 'brown'), c('Liver', 'cyan'), c('Liver', 'tan'))) {

	# get relevant strings to access each object
	tissue = thing[[1]]
	module = thing[[2]]
	tcode = substr(tissue, 1, 1)

	# get all the data we want in this table
	# start with module membership
	out.table = get(paste0(tolower(tcode), '.geneModuleMembership')) %>% 
				# pick desired module 
				select(paste0('MM_', tcode, '_', module)) %>%
				# get gene ids as own column
				mutate(gene_id = rownames(get(paste0(tolower(tcode), '.geneModuleMembership')))) %>%
				# add on gene names
				left_join(annot %>% select(gene_id, gene_name), 
							by = 'gene_id') %>%
				# add on module assignments
				left_join(get(paste0(tolower(tcode), '.moduleassign')), by = 'gene_id') %>%
				# take only the genes in this module
				filter(get(paste0(tolower(tissue), '_module')) == module) %>%
				# reorder columns nicely
				select(gene_id, gene_name, paste0('MM_', tcode, '_', module)) %>%
				# crop all numbers to 3 decimal places
				mutate(across(where(is.numeric), round, 3)) %>%
				# arrange by module membership
				arrange(get(paste0('MM_', tcode, '_', module))) %>%
				# give the columns sleek names
				rename(`Gene` = gene_name,
						`Ensembl ID` = gene_id,
						!!paste(tissue, module, 'membership (r)') := paste0('MM_', tcode, '_', module))
	# output to a csv excel will like
	write.csv(out.table, paste0('../../data/supplement_table_wgcna_', 
								tolower(tissue), '_', module, '.csv'),
								row.names = F)
}




