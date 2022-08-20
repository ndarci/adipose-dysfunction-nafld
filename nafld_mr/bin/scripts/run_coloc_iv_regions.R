## runColoc.R
# runs coloc (Bayesian form) to determine colocalization of 
# cis-eQTLs and TG/nonNAFLD GWAS variants in IV regions
library(coloc)

# import FPKM data, cropped to just the candidate region genes
fpkmpath = '../../data/KOBS_Invrs_FPKM_Baseline_All_peer25_invnorm_IVregions.csv'
fpkm = read.table(fpkmpath, header = TRUE, row.names = 1, sep = ',')

# set constants
gwasthresh = 5e-8
gwasN = 9926107

ppthresh = 0.8

# for each serum biomarker candidate, run coloc
# keep track of significant snps in iv candidate regions
runColoc <- function(gene) {
    # import eqtl data
    eqtlpath = paste('../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Adipose_cis_', gene, '.txt', sep = '')
    eqtl = read.table(eqtlpath, header = TRUE)

    # import overlapped eQTL-GWAS data
    ov = read.table(paste0('../../data/snp_overlap/TG_noNAFLD_', gene, '_SNPoverlap.txt'), header = T)

    # calculate var(beta)
    ov['varbeta_eqtl'] = (ov$beta_eqtl / ov$tstat_eqtl)^2
    ov['varbeta_gwas'] = (ov$se_gwas)^2

    # calculate sd(FPKM)
    sdfpkm = sd(fpkm[rownames(fpkm) == gene,])

    # run coloc
    # dataset1 = eQTL data
    # dataset2 = GWAS data
    resultTG <- coloc.abf(dataset1 = list(type = 'quant',
                                sdY = sdfpkm,
                                beta = ov$beta_eqtl,
                                varbeta = ov$varbeta_eqtl,
                                N = nrow(eqtl)),
                        dataset2 = list(type = 'quant',
                                beta = ov$beta_gwas,
                                varbeta = ov$varbeta_gwas,
                                N = gwasN),
                        MAF = ov$a1freq_gwas)

    ov.sig = NA
    # if significant colocalization, add significant overlapping SNPs to IV table
    if (resultTG$summary['PP.H4.abf'] > ppthresh) {
        # select significant SNPs to be used as IVs
        ov.sig = ov[ov$fdr_eqtl < 0.05 | ov$p_gwas < gwasthresh,]
    }

    return(list(res = resultTG, snps = ov.sig))
}

cols = c('gene_ID', 'nsnps', 'PP.H0.abf', 'PP.H1.abf',
            'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf')
outTG = data.frame(matrix(nrow = nrow(fpkm), ncol = length(cols)))
names(outTG) = cols

ivsnps = data.frame()

# iterate over all candidate region genes
for (i in 1:length(rownames(fpkm))) {
    id = rownames(fpkm)[i]

    # test for colocalization
    colocresult <- runColoc(id)

    # extract both output objects
    result = colocresult$res
    snps = colocresult$snps

    if (!is.na(snps)) {
        # add to overall dataframe
        ivsnps = rbind(ivsnps, snps)
    }

    # concat to result dataframe
    outTG[i,] = c(id, result$summary)
}

# convert numeric columns to numeric
num.cols = c('nsnps', 'PP.H0.abf', 'PP.H1.abf',
        'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf')
outTG[num.cols] <- sapply(outTG[num.cols], as.numeric)

# add gene symbol
# read in and clean up gene name and location information
annot = read.table('/u/project/pajukant/nikodm/kobs_limmaVoom/data/gencodeV19_annotations_formatted_genomewide.txt', header = TRUE, sep = '\t')
colnames(annot) = c('gene_ID', 'gene_ID_GENCODEver', 'chromosome', 'feature', 'startPos', 'endPos', 'strand', 'gene_type', 'gene_symbol')
annot = annot[,c('gene_ID', 'gene_symbol')]

outTG = merge(outTG, annot, by = 'gene_ID', all.x = T)

write.table(outTG, file = '../../data/coloc_IVregion_cis-eQTL_UKBgwasTG.txt')

# select significant colocalizations
outTG.sig = outTG[outTG$PP.H4.abf > ppthresh,]
write.table(outTG.sig, file = paste0('../../data/coloc_IVregion_cis-eQTL_UKBgwasTG_', ppthresh, 'PPH4.txt'))

# write list of candidate IV snps
write.table(ivsnps, '../../data/iv_snps_colocalized.txt')





