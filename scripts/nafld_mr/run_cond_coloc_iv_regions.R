# runs conditional coloc (Bayesian form) to determine colocalization of 
# cis-eQTLs and TG/nonNAFLD GWAS variants in IV regions
library(coloc)
library(tidyverse)

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

direc = commandArgs(trailingOnly = T)[1]

# read in relevant data
if (direc == 'f') {
    fpkmpath = '../../data/KOBS_Invrs_FPKM_Baseline_All_peer25_invnorm_IVregions.csv'
    fpkm = read.table(fpkmpath, header = TRUE, row.names = 1, sep = ',')
    eqtlprefix = '../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Adipose_cis_'
    ovtag = 'TG_noNAFLD_'
} else if (direc == 'b') {
    fpkmpath = '../../data/KOBS_liver_hg19_FPKM_nonzero_INT_correctedID_IVregions_backward.csv'
    fpkm = read.table(fpkmpath, header = TRUE, row.names = 1, sep = ',')
    rownames(fpkm) = stripVersion(rownames(fpkm))
    eqtlprefix = '../../data/ciseqtl_gene/KOBS_eQTL_BaselineAll_Liver_cis_'
    ovtag = 'NAFLD_noTG_'
} else {
    print("Invalid input argument")
}

# set constants
gwasthresh = 5e-8
gwasN = 9926107

ppthresh = 0.8

# for each serum biomarker candidate, run coloc
# keep track of significant snps in iv candidate regions
runCondColoc <- function(gene) {
    print(gene)
    # import eqtl data
    eqtlpath = paste(eqtlprefix, gene, '.txt', sep = '')
    eqtl = read.table(eqtlpath, header = TRUE)

    # import overlapped eQTL-GWAS data
    ov = read.table(paste0('../../data/snp_overlap/', ovtag, gene, '_SNPoverlap.txt'), header = T)

    # calculate var(beta)
    ov['varbeta_eqtl'] = (ov$beta_eqtl / ov$tstat_eqtl)^2
    ov['varbeta_gwas'] = (ov$se_gwas)^2

    # calculate sd(FPKM)
    sdfpkm = sd(fpkm[rownames(fpkm) == gene,])

    # import regional LD matrix
    ld = read.table(paste0('../../data/ld_matrix/ld_matrix_', gene, '.ld'), header = F)

    # import SNPs in the LD matrix (overlap of KOBS, UKB, METSIM)
    snplist = read.table(paste0('../../data/ld_matrix/ld_matrix_', gene, '.snplist'), header = F)$V1
    rssplit = strsplit(snplist, ':', fixed = T)
    snplist_rsonly = sapply(rssplit, '[[', 1)

    # give the ld matrix its row/col names
    rownames(ld) = snplist_rsonly
    names(ld) = snplist_rsonly

    # crop overlap table and LD matrix to a common list of SNPs
    snplist_rsonly = intersect(snplist_rsonly, ov$rsID_gwas)
    ov = ov %>% filter(rsID_gwas %in% snplist_rsonly)
    ld = ld[snplist_rsonly, snplist_rsonly]

    # verify that we're working with the exact same list in both tables
    all(ov$rsID_gwas == snplist_rsonly)

    # convert to matrix, not dataframe
    ld = as.matrix(ld)

    # prep for coloc run
    # dataset1 = eQTL data
    d1_eqtl = list(type = 'quant',
                    sdY = sdfpkm,
                    beta = ov$beta_eqtl,
                    varbeta = ov$varbeta_eqtl,
                    N = nrow(eqtl),
                    snp = ov$rsID_gwas)

    # dataset2 = GWAS data
    d2_gwas = list(type = 'quant',
                    beta = ov$beta_gwas,
                    varbeta = ov$varbeta_gwas,
                    N = gwasN,
                    snp = ov$rsID_gwas)

    # run conditional coloc
    result <- coloc.signals(method = 'cond',
                        mode = 'iterative',
                        LD = ld,
                        p12 = 1e-5,
                        MAF = ov$a1freq_gwas,
                        dataset1 = d1_eqtl,
                        dataset2 = d2_gwas)

    # check for significant colocalization between pairs of SNPs
    res.sig = result$summary %>% filter(PP.H4.abf > ppthresh)

    # # if significant colocalization, add significant overlapping SNPs to IV table
    # ov.sig = NA
    # if (nrow(res.sig) > 0) {
    #     # select significant SNPs to be used as IVs
    #     ov.sig = ov %>% filter(fdr_eqtl < 0.05 | p_gwas < gwasthresh)
    # }

    # return(list(res = res.sig, snps = ov.sig))
    return(list(res = res.sig))
}

cols = c('gene_ID', 'hit2', 'hit1', 'nsnps', 'PP.H0.abf', 'PP.H1.abf',
            'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf')
outTG = data.frame(matrix(nrow = 0, ncol = length(cols)))
names(outTG) = cols

# ivsnps = data.frame()

# iterate over all candidate region genes
for (id in rownames(fpkm)) {
    # test for colocalization
    colocresult <- runCondColoc(id)

    # extract both output objects
    result = colocresult$res %>% data.frame(.)
    # snps = colocresult$snps

    # if there was colocalization
    if (nrow(result > 0)) {
        # add to overall dataframes
        # ivsnps = rbind(ivsnps, snps)

        result['gene_ID'] = id

        # concat to result dataframe
        outTG = rbind(outTG, result %>% select(gene_ID, hit2, hit1, nsnps, PP.H0.abf, PP.H1.abf, 
                                                    PP.H2.abf, PP.H3.abf, PP.H4.abf))
    }
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

# select significant colocalizations
# print(outTG)
write.table(outTG, file = paste0('../../data/cond_coloc_IVregion_cis-eQTL_UKBgwasTG_', 
                                    ppthresh, 
                                    'PPH4',
                                    rep('_backward', direc == 'b'),
                                    '.txt'))

# # write list of candidate IV snps
# write.table(ivsnps, '../../data/iv_snps_cond_colocalized.txt')







