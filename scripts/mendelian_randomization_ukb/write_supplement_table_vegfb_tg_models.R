library(tidyverse)

# read in model performance stats
log = read.table('../../data/vegfb_mr_snp_tg_models_logistic.txt', header = T)
lin = read.table('../../data/vegfb_mr_snp_tg_models_linear.txt', header = T)
ela = read.table('../../data/vegfb_mr_snp_tg_models_elastic.txt', header = T)

# global vectors the script will use
keepvars = c('Triglycerides', 'VEGFB_INT+Triglycerides')
hist_order = c('steatosis', 'fibrosis', 'diagnosis')

clean_up_table <- function(rawdf) {
	cleandf = rawdf %>%
		# just keep the models we care about
		filter(vars_in_model %in% keepvars) %>%
		# fix the order of models in vars in model and phenotype
		mutate(vars_in_model = factor(vars_in_model, levels = keepvars),
				phenotype = factor(phenotype, levels = hist_order)) %>%
		arrange(phenotype, vars_in_model) %>%
				# take '_INT' off VEGFB
		mutate(vars_in_model = recode(vars_in_model,
				`VEGFB_INT+Triglycerides` = 'VEGFB+Triglycerides'),
				# make the phenotypes look pretty
				phenotype = recode(phenotype,
						steatosis = 'Steatosis',
						fibrosis = 'Fibrosis',
						diagnosis = 'NASH')) %>%
		# convert numeric to 3 decimal points
		mutate(across(where(is.numeric), round, 3))
	return(cleandf)
}

# clean up logistic table
log.clean = clean_up_table(log) %>% 
	# make the model type capitalized
	mutate(model_type = recode(model_type,
			logistic = 'Logistic')) %>%
	# make the col names nice
	rename(Phenotype = phenotype,
			`Variables in model` = vars_in_model,
			`Model type` = model_type,
			`Pseudo-r2` = pseudo.r2.Nagelkerke)

# clean up linear table
lin.clean = clean_up_table(lin) %>%
	# make the model type capitalized
	mutate(model_type = recode(model_type,
			linear = 'Linear')) %>%
	# make the col names nice
	rename(Phenotype = phenotype,
		`Variables in model` = vars_in_model,
		`Adj. r2` = adj.r2,
		`Model type` = model_type)

# clean up elastic table
ela.clean = clean_up_table(ela) %>%
	# remove SNP column
	select(-c(rs2845885_coef)) %>%
	# make the model type capitalized
	mutate(model_type = recode(model_type,
			elasticnet = 'Elastic Net')) %>%
	# make the col names nice
	rename(Phenotype = phenotype,
		`Variables in model` = vars_in_model,
		`VEGFB coeff.` = VEGFB_coef,
		`TG coeff.` = TG_coef,		
		`Model type` = model_type)

write.table(log.clean, '../../data/supplement_table_vegfb_tg_logistic.txt', row.names = F)
write.table(lin.clean, '../../data/supplement_table_vegfb_tg_linear.txt', row.names = F)
write.table(ela.clean, '../../data/supplement_table_vegfb_tg_elastic.txt', row.names = F)




