library(tidyverse)

# read in DE results
amt_order = c('0ng', '20ng')
# prot_order = c('CCDC80')
prot_order = c('CCDC80', 'SOD3')
read_in_de <- function(nominal = F) {
	de = data.frame()
	for(protname in prot_order) {
		covcorrect = F 
		covariate = ""
		fn = paste0("../../data/topTable_", 
					protname, 
					rep(paste0("_covcorrected_", covariate), covcorrect), 
					rep("_nominal", nominal),
					".txt")
		newdata = read.table(fn, header = T, sep = '\t')
		if (nrow(newdata) > 0) {
			newdata["protein_name"] = protname
			newdata["covcorrect"] = covcorrect
			newdata["covariate"] = covariate
			de = rbind(de, newdata)
		}
	}
	return(de)
}

de = read_in_de()

# read in Ilakya's DESeq2 results to compare against mine
deseq = read.table('../../data/CCDC80_RNAseq_result_for_Niko/DESeq2_result/2023-03-24_20-50-05_featureCounts-20ng_VS_0ng-DESeq2-result.txt', header = T)

# convert EU format numbers to USA format
eu_to_usa = function(eu_number) {
	usa_number = as.numeric(gsub(',', '.', eu_number))
	return(usa_number)
}
deseq = deseq %>%
	mutate(across(c(baseMean, log2FoldChange, 
		lfcSE, stat, pvalue, padj), eu_to_usa)) %>%
	# select only significant genes
	filter(padj < 0.05) %>%
	# drop cols we don't care about
	select(geneName, log2FoldChange, pvalue, padj) %>%
	rename('gene_name' = geneName,
			'logFC_ilakya' = log2FoldChange,
			'pval_ilakya' = pvalue,
			'padj_ilakya' = padj)

# check for overlap between my DE genes and Ilakya's
ov = de %>% 
	select(gene_ID, gene_name, protein_name, logFC, P.Value, adj.P.Val) %>%
	inner_join(deseq, by = 'gene_name')
print("overlap with Ilakya's DE genes:")
print(ov)


# read in zong's NAFLD_imp GWAS results to see if any of the DE genes have nearby variants associated with NAFLD
gwas = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen_5e-8.txt', header = T)


# read in gene position data
annot = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_formatted_genomewide.txt', header = TRUE, sep = '\t')
colnames(annot) = c('gene_ID', 'gene_ID_GENCODEver', 'chromosome', 'feature', 'startPos', 'endPos', 'strand', 'gene_type', 'gene_symbol')
annot = annot %>%
        select(gene_symbol, chromosome, startPos, endPos) %>%
        mutate(chr = substr(chromosome, 4, nchar(chromosome))) %>%
        rename(gene_name = gene_symbol)

# merge position onto DE data
de = de %>% left_join(annot, by = 'gene_name')

# find all gwas hits in each DE gene's cis region
# define the cis region
cisdist = 5e5
de = de %>% mutate(cis_start = startPos - cisdist,
					cis_end = endPos + cisdist)
de[de$cis_start < 0, 'cis_start'] <- 0 # bump ones that go below 0 back to 0


de.gwas = data.frame()
for (gene in de$gene_name) {
	# get the region we're searching in on the same chr as the gene
	this_chr = de[de$gene_name == gene, 'chr']
	this_cis_start = de[de$gene_name == gene, 'cis_start']
	this_cis_end = de[de$gene_name == gene, 'cis_end']

	# select all significant GWAS hits in this region
	gwas.sub = gwas %>% filter(CHR == this_chr && 
								BP > this_cis_start && 
								BP < this_cis_end) %>%
						select(SNP, CHR, BP, A1FREQ, BETA, P_BOLT_LMM_INF) %>%
						mutate(gene_region = gene)

	# if there were any hits, add it to the overall table
	if(nrow(gwas.sub) > 0) {
		de.gwas = rbind(de.gwas, gwas.sub)
	}
}
print("NAFLD gwas hits in the cis regions of DE genes")
print(de.gwas)





