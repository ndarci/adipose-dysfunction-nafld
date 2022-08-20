library(tidyverse)

stripVersion <- function(idlist) {
	strid = strsplit(idlist, ".", fixed = T)
	stripped = sapply(strid, "[[", 1)
	return(stripped)
}

# read in correlation results
c = read.table('../../data/sbc_all_liver_gene_correlations.txt', header = T)

# # overlap with wnt pathway genes
# wnt = read.table('/u/project/pajukant/nikodm/sbc_sirna/data/wnt_pathway_genes.txt')$V1

# wnt_sig = c[c$othergene_liver_sym %in% wnt & c$pvalue < 0.05,] %>% na.omit(.)

# summarise cross-tissue correlation results
bonfthresh = 0.05/nrow(c)
c.sum = c %>% group_by(sbc_adipose_sym) %>%
	filter(pvalue < bonfthresh) %>%
	summarise(n_cor = length(othergene_liver),
			n_pos_cor = length(othergene_liver[r > 0]),
			n_neg_cor = length(othergene_liver[r < 0]))

write.table(c.sum, '../../data/xtissue_ligand_cor_per_sbc.txt', row.names = F, quote = F)

# export SFRP2 correlated genes
sfrp2cor = c %>% filter(sbc_adipose_sym == 'SFRP2' & pvalue < bonfthresh) %>% 
			.$othergene_liver
# sfrp2cor_sym = c %>% filter(sbc_adipose_sym == 'SFRP2' & pvalue < bonfthresh) %>% 
# 			.$othergene_liver_sym
write.table(sfrp2cor, '../../data/SFRP2_correlated_liver_genes.txt', col.names = F, row.names = F, quote = F)
# write.table(sfrp2cor_sym, '../../data/SFRP2_correlated_liver_genes_symbol.txt', col.names = F, row.names = F, quote = F)

# export all significant correlations
c.sig = c %>% 
	filter(pvalue < bonfthresh)

write.table(c.sig %>% select(othergene_liver) %>% unique, '../../data/allSBC_correlated_liver_genes.txt', col.names = F, row.names = F, quote = F)

# check if correlated genes are DE in liver
# read in liver DE results
liverde = data.frame()
for (pheno in c('steatosis', 'fibrosis', 'diagnosis')) {
	newdata = read.table(paste0('/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/topTable_', pheno, '_liver.txt'), header = T)
	newdata['gene_ID'] = stripVersion(rownames(newdata))
	newdata['pheno'] = pheno
	liverde = rbind(liverde, newdata)
}

# group on gene ID and filter for significance
# pivot to get all 3 phenotype columns separately
liverde = liverde %>% 
			filter(adj.P.Val < 0.05) %>%
			group_by(gene_ID) %>%
			summarise(pheno = pheno,
						logFC = logFC,
						adjP = adj.P.Val) %>%
			# mutate(foo = T) %>%
			pivot_wider(names_from = pheno, values_from = c(logFC, adjP))
names(liverde)[2:ncol(liverde)] = paste0('liverDE_', names(liverde)[2:ncol(liverde)])

# merge with correlation results
c.sig = merge(c.sig, liverde, by.x = 'othergene_liver', by.y = 'gene_ID', all.x = T)

# check which WGCNA module each of the correlated genes belongs to
# import gene-module matching info
lnames = load('../../data/xtissue_correlation_results.RData')

# merge onto existing table
c.sig = merge(c.sig, l.moduleassign, by.x = 'othergene_liver', by.y = 'gene_id') %>%
		rename(othergene_liver_module = 'liver_module') %>%
		merge(., a.moduleassign, by.x = 'sbc_adipose', by.y = 'gene_id') %>%
		rename(sbc_adipose_module = 'adipose_module')

# check how many GWAS hits are in each correlated gene's cis region  
# import gwas summary stats
gwas.nafld = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen_5e-2.txt', header = T)
gwas.alt = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/ALT/ALT_bgen.5e-2.txt', header = T)
# gwas.ast = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/ALT/ALT_bgen.5e-2.txt', header = T)
gwas.ggt = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/GGT/GGT_bgen.5e-2.txt', header = T)

# import annotations with gene coordinates
annot = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_formatted_genomewide.txt', header = TRUE, sep = '\t')
colnames(annot) = c('gene_ID', 'gene_ID_GENCODEver', 'chromosome', 'feature', 'startPos', 'endPos', 'strand', 'gene_type', 'gene_symbol')
annot = annot %>% 
	select(gene_ID, gene_symbol, chromosome, startPos, endPos) %>% 
	mutate(chr = substr(chromosome, 4, nchar(chromosome))) %>%
	filter(gene_ID %in% c.sig$othergene_liver)

gwasthresh = 5e-8
cisdist = 1e6

# iterate over all GWASes
for (gwastrait in c('nafld', 'alt', 'ggt')) {
	# make a column for number of cis GWAS hits at each gene
	annot[paste0('cis_gwas_hits_', gwastrait)] = NA
	# iterate over all correlated genes
	for (i in 1:nrow(annot)) {
		# point to the gwas results for this trait
		gwas = get(paste0('gwas.', gwastrait))
		# count significant GWAS hits in this gene's cis region
		nhits = gwas %>% filter(CHR == annot[i,'chr'] &
								BP > (annot[i,'startPos']-cisdist) &
								BP < (annot[i,'endPos']+cisdist) &
								P_BOLT_LMM_INF < gwasthresh) %>%
						nrow(.)
		# add result to gwas count column
		annot[i, paste0('cis_gwas_hits_', gwastrait)] <- nhits
	}
}

# merge back onto original table
c.sig = merge(c.sig, annot %>% select(gene_ID, starts_with('cis_gwas_hits')),
				by.x = 'othergene_liver', by.y = 'gene_ID', all.x = T)

# add info on whether the correlation was nominally significant in healthy individuals
# read in healthy-only correlation results
ch = read.table('../../data/sbc_all_liver_gene_correlations_healthy.txt', header = T)

# summarise healthy-only correlation results
ch.sum = ch %>% group_by(sbc_adipose_sym) %>%
	filter(pvalue < 0.05) %>%
	summarise(n_cor = length(othergene_liver),
			n_pos_cor = length(othergene_liver[r > 0]),
			n_neg_cor = length(othergene_liver[r < 0]))

ch.sig = ch %>% filter(pvalue < 0.05) %>%
		select(-ends_with('_sym'))

names(ch.sig) = paste0('healthy_', names(ch.sig))
c.sig = merge(c.sig, ch.sig, 
			by.x = c('sbc_adipose', 'othergene_liver'),
			by.y = c('healthy_sbc_adipose', 'healthy_othergene_liver'),
			all.x = T)

c.sig['nom_sig_in_healthy'] = !is.na(c.sig$healthy_r)
table(c.sig$nom_sig_in_healthy)

c.sig = c.sig %>% select(
						sbc_adipose, othergene_liver, 
						sbc_adipose_sym, othergene_liver_sym,
						r, pvalue, r2,
						starts_with('liverDE_'),
						sbc_adipose_module, othergene_liver_module,
						starts_with('cis_gwas_hits_'),
						nom_sig_in_healthy,
						starts_with('healthy_')
						)

write.table(c.sig, '../../data/xtissue_ligand_cor_stats.txt', row.names = F, quote = F)




