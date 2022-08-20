library(tidyverse)

direc = commandArgs(trailingOnly = T)[1]

# read in and clean up gene name and location information
annot = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_formatted_genomewide.txt', header = TRUE, sep = '\t')
colnames(annot) = c('gene_ID', 'gene_ID_GENCODEver', 'chromosome', 'feature', 'startPos', 'endPos', 'strand', 'gene_type', 'gene_symbol')
annot = annot %>% 
	select(gene_ID, gene_symbol, chromosome, startPos, endPos) %>% 
	mutate(chr = substr(chromosome, 4, nchar(chromosome)))

# read in relevant data
if (direc == 'f') {
	# adipose aware DE gene list
	de = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/DE_genes_filter_except_secreted_all.txt')
	# exposure GWAS (TG)
	gwas.exp = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/TG/TG_bgen.5e-2.txt', header = T)
	# outcome GWAS (NAFLD)
	gwas.out = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen_5e-2.txt', header = T)
	# set eqtl file location prefix for reading in eqtl data later
	eqtlprefix = '../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Adipose_cis_'
	# set important part of output file name
	outtag = 'TG_noNAFLD_'
} else if (direc == 'b') {
	# liver aware DE gene list
	de = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/DE_genes_liver_not_adipose_all.txt')
	# exposure GWAS (NAFLD)
	gwas.exp = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/NAFLD_imp/NAFLD_imp_bgen_5e-2.txt', header = T)
	# outcome GWAS (TG)
	gwas.out = read.table('/u/project/pajukant/zmiao/UKBiobank/GWAS/TG/TG_bgen.5e-2.txt', header = T)
	# eqtl file prefix for liver eqtls
	eqtlprefix = '../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Liver_cis_'
	# set important part of output file name
	outtag = 'NAFLD_noTG_'
} else {
	print("Invalid input argument")
}
names(de)[1] = 'gene_ID'

# add gene positions and define cis regions
cisdist = 1e6
de = de %>% 
		unique(.) %>% 
		merge(., annot, by = 'gene_ID', all.x = T) %>%
		mutate(cis_start = startPos - cisdist,
				cis_end = endPos + cisdist)
de[de$cis_start < 0, 'cis_start'] <- 0 # fix start positions that dip below 0

# remove outcome gwas variants from consideration
gwas_thresh = 5e-8
gwas.exp = gwas.exp %>% mutate(chrbp = paste0(CHR, ':', BP))
gwas.out.sig = gwas.out %>% mutate(chrbp = paste0(CHR, ':', BP)) %>%
							filter(P_BOLT_LMM_INF < gwas_thresh)
gwas.exp = gwas.exp %>% filter(!(chrbp %in% gwas.out.sig$chrbp))

# add a column for the IV selection
de['iv'] = NA
de['n_snps'] = NA

# scan across all input genes to check if they have a cis-eQTL and GWAS hit in their cis region
for (thisgene in de$gene_ID) {
	# thisgene = 'ENSG00000008226'

	# import gene-specific eqtl data 
	eqtl = read.table(paste0(eqtlprefix, thisgene,'.txt'), header = T)
	if (nrow(eqtl) == 0) {
		de[de$gene_ID == thisgene, 'iv'] <- F 
		next
	}
	strsnp = strsplit(eqtl$SNP, ':', fixed = T)

	# add chromosome and genetic position columns 
	# different for different snp ID formats
	if (direc == 'f') {
		eqtl['chr'] = sapply(strsnp, '[[', 1)
		eqtl['pos'] = sapply(strsnp, '[[', 2)
	} else if (direc == 'b') {
		eqtl['rsID'] = sapply(strsnp, '[[', 1)
		eqtl = eqtl[grepl('rs', eqtl$rsID),]
		eqtl = merge(eqtl, gwas.exp %>% 
							select(SNP, CHR, BP) %>% 
							rename(rsID = 'SNP')) %>%
					rename(chr = 'CHR', pos = 'BP')
	} else {
		print("Invalid input argument")
	}

	# restrict eqtl and gwas to cis distance from this gene
	thischrom = de[de$gene_ID == thisgene, 'chr']
	thiscis_start = de[de$gene_ID == thisgene, 'cis_start']
	thiscis_end = de[de$gene_ID == thisgene, 'cis_end']

	eqtl.sub = eqtl %>% filter(chr == thischrom &
								pos > thiscis_start & 
								pos < thiscis_end)
	eqtl.sub.sig = eqtl.sub %>% filter(FDR < 0.05)

	gwas.sub = gwas.exp %>% filter(CHR == thischrom &
								BP > thiscis_start &
								BP < thiscis_end)
	gwas.sub.sig = gwas.sub %>% filter(P_BOLT_LMM_INF < gwas_thresh)

	# check if this gene goes into the IV list
	yesiv = nrow(eqtl.sub.sig) > 0 & nrow(gwas.sub.sig) > 0

	# update the table accordingly
	de[de$gene_ID == thisgene, 'iv'] <- yesiv

	# if it's an IV region, write ALL overlapping SNPs in the region to a new table
	if(yesiv == T) {
		# format eqtl table with desired features
		eqtl.sub = eqtl.sub %>% 
			select(gene, beta, t.stat, p.value, FDR, chr, pos) %>%
			rename(gene_eqtl = gene,
					beta_eqtl = beta,
					tstat_eqtl = t.stat,
					p_eqtl = p.value,
					fdr_eqtl = FDR,
					chr_eqtl = chr,
					pos_eqtl = pos)

		# format gwas table with desired features
		gwas.sub = gwas.sub %>%
			select(BETA, SNP, ALLELE1, ALLELE0, A1FREQ, SE, P_BOLT_LMM_INF, CHR, BP) %>%
			rename(beta_gwas = BETA,
					rsID_gwas = SNP,
					a1_gwas = ALLELE1,
					a0_gwas = ALLELE0,
					a1freq_gwas = A1FREQ,
					se_gwas = SE,
					p_gwas = P_BOLT_LMM_INF,
					chr_gwas = CHR,
					pos_gwas = BP)

		# merge on snp position
		ov = merge(eqtl.sub, gwas.sub,
					by.x = c('chr_eqtl', 'pos_eqtl'),
					by.y = c('chr_gwas', 'pos_gwas'))

		# remove variants with no rsID
		ov = ov[grepl('rs', ov$rsID_gwas),]

		# remove multiallelic variants
		ov = ov[!(duplicated(ov$rsID_gwas) | duplicated(ov$rsID_gwas, fromLast = T)),]

		# remember number of SNPs in this region
		de[de$gene_ID == thisgene, 'n_snps'] <- nrow(ov)

		# save the overlap file
		write.table(ov, paste0('../../data/snp_overlap/', outtag, thisgene, '_SNPoverlap.txt'))
		
		# generate a simpler table with only the SNP positions for LD matrix
		snplist = ov %>% select(rsID_gwas, a1_gwas, a0_gwas)

		write.table(snplist, 
					paste0('../../data/ld_matrix/', outtag, thisgene, '_SNPoverlap_SNPlist.txt'), 
					sep = ':', row.names = F, col.names = F, quote = F)
	}
}

# output a list of IV SNP regions
table(de$iv)
de.iv = de %>% filter(iv == T)
print(de.iv)

write.table(de.iv, paste0('../../data/iv_candidate_regions', 
							rep('_backward', direc == 'b'),
							'.txt'))
write.table(de.iv$gene_ID, paste0('../../data/iv_candidate_genes_ensembl',
							rep('_backward', direc == 'b'),
							'.txt'), 
							quote = F, row.names = F, col.names = F)







