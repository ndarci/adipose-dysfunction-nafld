library(tidyverse)

# read in correlation results
c = read.table('../../data/sbc_all_liver_module_correlations.txt', header = T)

# get bonferroni corrected significance
c = c %>% mutate(pvalue_adj = pvalue * nrow(c),
				pvalue_adj = ifelse(pvalue_adj > 1, 1, pvalue_adj))

# select only relevant cols
c = c %>% select(sbc_adipose_sym, othergene_liver, r, r2, pvalue, pvalue_adj) %>%
		rename(liver_module = othergene_liver)

# get significant correlations
c.sig = c %>% 
		filter(pvalue < 0.05) %>%
		arrange(sbc_adipose_sym, r)

# get bonf significant correlations
c.sig.bonf = c.sig %>% filter(pvalue_adj < 0.05)

write.table(c.sig, '../../data/adipose_sbc_liver_module_corr_nominalsig.txt', row.names = F)
write.table(c.sig.bonf, '../../data/adipose_sbc_liver_module_corr_bonferronisig.txt', row.names = F)
