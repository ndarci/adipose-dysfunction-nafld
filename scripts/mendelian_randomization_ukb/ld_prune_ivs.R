library(tidyverse)

# import ld between ivs and nafld variants
ld = read.table('../../data/ld_matrix_ivs_nafldhits.ld')
# collect snp names
snplist = read.table('../../data/ld_matrix_ivs_nafldhits.snplist')$V1
rssplit = strsplit(snplist, ':', fixed = T)
snplist_rsonly = sapply(rssplit, '[[', 1)
names(ld) = snplist_rsonly
rownames(ld) = snplist_rsonly

# import current list of IVs
iv = read.table('../../data/iv_snps_cond_colocalized_clean.txt')

# get list of nafld variants
nafld = names(ld)[!(names(ld) %in% iv$hit1)]

# clip ld matrix to just ivs and nafld variants
iv_nafld = ld[iv$hit1, nafld]

# check if any correlations have r2>0.2
ld_thresh = 0.2
print('correlations with NAFLD GWAS hits:')
which(iv_nafld > ld_thresh, arr.ind = T) %>% 
					data.frame %>% 
					arrange(row)

# clip ld matrix to just IVs
iv_iv = ld[iv$hit1, iv$hit1]

# plot heatmap of correlations
# get data into correct format
cor.melt = iv_iv %>% 
			mutate(iv1 = rownames(iv_iv)) %>% 
			pivot_longer(cols = starts_with('rs'), names_to = 'iv2', values_to = 'r2')

# generate the plot
ivcorr = ggplot(cor.melt, aes(x = iv1, y = iv2, fill = r2)) +
    geom_tile() + 
    geom_text(aes(label = round(r2, 2))) +
    coord_fixed() +
    scale_fill_gradient(low = 'white', 
                        high = 'orangered') +
    # guides(fill = guide_colorbar(title = 'r2',
    #                             barwidth = 0.5, 
    #                             barheight = 15,
    #                             label.position = 'right',
    #                             title.hjust = 0,
    #                             ticks = F)) + 
    # scale_y_discrete(limits = pc1_ord, position = 'left') +
    # scale_x_discrete(limits = pc1_ord) +
    theme(panel.background = element_blank(),
            legend.position = 'right',
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 14))

ggsave('../../figs/ldmatrix_iv_candidates.png', ivcorr)

# identify which IVs are too correlated
iv.ld = which(iv_iv > ld_thresh, arr.ind = T) %>% 
					data.frame %>% 
					filter(row != col) %>%
					arrange(row) %>%
                    mutate(hit1 = rownames(iv_iv)[row])

# select representative IV of those in LD using cis and GWAS stats 
iv.stat = read.table('../../data/iv_snps_cond_colocalized_cis_gwas_stats.txt')

# import cis and gwas stats
iv.stat.ld = iv.stat %>% group_by(hit1) %>%
                        summarise(gene_ID = gene_ID,
                                fdr_eqtl = fdr_eqtl,
                                p_gwas = p_gwas) %>%
                        filter(hit1 %in% iv.ld$hit1)

# find snp with minimum p-values in each one
sig.cis = iv.stat.ld[which.min(iv.stat.ld$fdr_eqtl),]$hit1
sig.gwas = iv.stat.ld[which.min(iv.stat.ld$p_gwas),]$hit1

if(sig.cis != sig.gwas) {
    print('single SNP does not have most signif cis FDR and GWAS pvalue')
} else {
    keepIV = sig.cis

    # select the set of IVs that pass LD checks and go to the MR analysis
    iv.out = iv %>% filter(!(hit1 %in% iv.ld$hit1) | hit1 == keepIV)
    write.table(iv.out, '../../data/iv_snps_cond_colocalized_ldpruned.txt')

    write.table(iv.out$hit1, '../../data/iv_snps_cond_colocalized_ldpruned_just_rsID.txt', row.names = F, col.names = F, quote = F)
}














