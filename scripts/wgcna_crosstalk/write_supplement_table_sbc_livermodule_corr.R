library(tidyverse)

# read in correlation data
c = read.table('../../data/adipose_sbc_liver_module_corr_bonferronisig.txt', header = T)

# clean everything up to make it nice for the paper
c.clean = c %>% 
		# get rid of L_ME in liver module names
		mutate(liver_module = substr(liver_module, 5, nchar(liver_module))) %>%
		# convert p-values to scientific notation
		mutate_at(c('pvalue', 'pvalue_adj'), formatC, format = 'e', digits = 3) %>%
		# crop all numbers to 3 decimal places
		mutate(across(where(is.numeric), round, 3)) %>%
		# sleek column names with no underscores
		rename(`SBC` = sbc_adipose_sym,
				`Liver module` = liver_module,
				`P-value` = pvalue,
				`Adj. p-value` = pvalue_adj)

write.csv(c.clean, '../../data/supplement_table_sbc_liver_module_corr.csv', row.names = F)