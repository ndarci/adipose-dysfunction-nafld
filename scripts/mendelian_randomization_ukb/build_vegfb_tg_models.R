library(tidyverse)
library(DescTools)
library(caret)
library(glmnet)
library(pROC)

read_clean_data <- function(pheno) {
	# import cleaned data from best subsets
	expr_cov = read.table(paste0('/u/project/pajukant/nikodm/kobs_limmaVoom/data/vegfb_logisticmodel_data_', pheno, '.txt'), header = T)
	expr_cov = expr_cov %>% rename(sample_ID = 'Row.names')

	# add TG to the data
	rawcov = read.csv('/u/project/pajukant/nikodm/kobs_limmaVoom/data/cov.csv')
	rawcov = rawcov %>% rename(sample_ID = 'X')
	df = expr_cov %>% left_join(rawcov %>% select(sample_ID, Triglycerides), by = 'sample_ID')
	# remove samples with TG missing
	df = na.omit(df)

	# add special MR VEGFB SNP genotypes to the data
	# the SNP is rs2845885
	geno = read.table('../../data/vegfb_mr_snp_genotypes_text.raw', header = T)
	geno = geno %>% 
		# make sample ID good for merging
		rename(sample_ID = 'IID') %>%
		mutate(sample_ID = paste0('X', sample_ID)) %>%
		# drop columns we don't care about
		select(sample_ID, rs2845885_C)
	df = df %>% 
		# merge onto the rest of the data
		left_join(geno, by = 'sample_ID') %>%
		# select the final set of columns
		select(sample_ID, paste0(pheno, '_grp'), VEGFB_INT, Triglycerides, rs2845885_C, 
			Age, isMale, RIN, UNIQ_MAP_PERCENT, PCT_INTRONIC_BASES, MEDIAN_3PRIME_BIAS) %>%
		# log10 transform the TG
		mutate(Triglycerides = log10(Triglycerides))

	return(df)
}

run_linear_model <- function(pheno, covariates, var_combo, df) {
	# set up the variables that will go in the model
	formula_string = paste(paste0(pheno, '_grp'), '~', 
				paste(c(var_combo, covariates), collapse = '+'))
	# fit predictors to the outcome 
	model = lm(as.formula(formula_string), data = df)
	# extract meaningful stats
	r2 = summary(model)$r.squared
	adj.r2 = summary(model)$adj.r.squared

	# pval = pf(summary(model)$fstatistic[1],
	# 		summary(model)$fstatistic[2],
	# 		summary(model)$fstatistic[3], lower.tail = F)

	# return a vector of stats from the model fit
	# return(c(r2, adj.r2, pval))
	return(c(r2, adj.r2))
}

run_logistic_model <- function(pheno, covariates, var_combo, df) {
	# set up the variables that will go in the model
	formula_string = paste(paste0(pheno, '_grp'), '~', 
				paste(c(var_combo, covariates), collapse = '+'))
	# fit predictors to the outcome 
	model = glm(as.formula(formula_string), data = df, family = binomial)
	# calculate pseudo r2
	pseudoR2 = PseudoR2(model, which = 'Nagelkerke')

	# calculate AUC
	auc = roc(df[,paste0(pheno, '_grp')], fitted(model))$auc

	# # calculate p-value using log likelihood ratio test
	# # get number of variables in model for degrees of freedom in xsq test
	# nvar_in_model = lengths(regmatches(formula_string, gregexpr("+", formula_string, fixed = T))) + 1
	# llr = -2 * (model$null.deviance - model$deviance)
	# pval = pchisq(llr, (nvar_in_model+1) - 1, lower.tail = FALSE)

	# return a vector of stats from the model fit
	# return(c(pseudoR2, pval))
	return(c(pseudoR2, auc))
}

run_elastic_model <- function(pheno, covariates, var_combo, df) {
	# elasticnet wants the outcome to be a factor or numeric
	df = df %>% 
		mutate("{paste0(pheno, '_grp')}" := as.factor(get(paste0(pheno, '_grp'))))

	# set the model formula
	formula_string = paste(paste0(pheno, '_grp'), '~', 
		paste(c(var_combo, covariates), collapse = '+'))

	# train the model to all the data using CV
	set.seed(22)
	elastic = train(as.formula(formula_string), data = df, 
			method = 'glmnet',
			trControl = trainControl('cv', number = 10),
			tuneLength = 10)

	# extract model coefficients
	elastic_coef = coef(elastic$finalModel, elastic$bestTune$lambda)
	vegfb_coef = ifelse('VEGFB_INT' %in% var_combo, elastic_coef[,1]['VEGFB_INT'], NA)
	tg_coef = ifelse('Triglycerides' %in% var_combo, elastic_coef[,1]['Triglycerides'], NA)
	snp_coef = ifelse('rs2845885_C' %in% var_combo, elastic_coef[,1]['rs2845885_C'], NA)

	# make predictions (using the training data/only data we have)
	predictions <- elastic %>% predict(df)

	# extract r2 of the fit
	Rsquare = R2(as.numeric(predictions), as.numeric(df[,paste0(pheno, '_grp')]))

	return(c(Rsquare, vegfb_coef, tg_coef, snp_coef))
}

# define the phenotypes we'll use as outcome variables
phenotypes = c('steatosis', 'fibrosis', 'diagnosis')

# define the different variable combos we'll compare in different models
# all models will correct for the RNA-seq covariates
var_combos = list(c('VEGFB_INT', 'Triglycerides', 'rs2845885_C'),
					c('VEGFB_INT', 'rs2845885_C'),
					c('VEGFB_INT', 'Triglycerides'),
					c('Triglycerides'),
					c('VEGFB_INT'),
					c('rs2845885_C'))

# set up the result table for the linear regression
# linear_stat_col_names = c('phenotype', 'vars_in_model', 'r2', 'adj.r2', 'pval')
linear_stat_col_names = c('phenotype', 'vars_in_model', 'r2', 'adj.r2')
linear_model_stats = data.frame(matrix(NA, nrow = length(phenotypes)*length(var_combos),
						ncol = length(linear_stat_col_names)))
names(linear_model_stats) = linear_stat_col_names

# set up the result table for the logistic regression
# logistic_stat_col_names = c('phenotype', 'vars_in_model', 'pseudo.r2', 'pval')
logistic_stat_col_names = c('phenotype', 'vars_in_model', 'pseudo.r2.Nagelkerke', 'AUC')
logistic_model_stats = data.frame(matrix(NA, nrow = length(phenotypes)*length(var_combos),
						ncol = length(logistic_stat_col_names)))
names(logistic_model_stats) = logistic_stat_col_names

# set up the same for the elasticnet regression
elastic_stat_col_names = c('phenotype', 'vars_in_model', 'r2', 'VEGFB_coef', 'TG_coef', 'rs2845885_coef')
elastic_model_stats = data.frame(matrix(NA, nrow = length(phenotypes)*length(var_combos),
						ncol = length(elastic_stat_col_names)))
names(elastic_model_stats) = elastic_stat_col_names

# i iterates over the rows in both result columns
i = 1
# covariates are the same in all models
covariates = c('Age', 'isMale', 'RIN', 'UNIQ_MAP_PERCENT', 
				'PCT_INTRONIC_BASES', 'MEDIAN_3PRIME_BIAS')
for (pheno in phenotypes) {
	# get the model data prepped for this phenotype
	df = read_clean_data(pheno)

	# build linear and logistic models comparing variance explained by VEGFB expression, TG, and rs2845885
	for (var_combo in var_combos) {
		# var_combo = var_combo[[1]]
		# fit linear model and update the linear model result table
		fit_stats = run_linear_model(pheno, covariates, var_combo, df)
		linear_model_stats[i,] = c(pheno, paste(var_combo, collapse = '+'), fit_stats)
		
		# do the same with the logistic model
		fit_stats_logistic = run_logistic_model(pheno, covariates, var_combo, df)
		logistic_model_stats[i,] = c(pheno, paste(var_combo, collapse = '+'), fit_stats_logistic)
		# print(i)
		# print(var_combo)
		# print(linear_model_stats)
		# print(logistic_model_stats)

		# add an elasticnet model as well
		fit_stats_elastic = run_elastic_model(pheno, covariates, var_combo, df)
		elastic_model_stats[i,] = c(pheno, paste(var_combo, collapse = '+'), fit_stats_elastic)

		i = i+1
	}
}
linear_model_stats['model_type'] = 'linear'
logistic_model_stats['model_type'] = 'logistic'
elastic_model_stats['model_type'] = 'elasticnet'

write.table(linear_model_stats, '../../data/vegfb_mr_snp_tg_models_linear.txt', row.names = F)
write.table(logistic_model_stats, '../../data/vegfb_mr_snp_tg_models_logistic.txt', row.names = F)
write.table(elastic_model_stats, '../../data/vegfb_mr_snp_tg_models_elastic.txt', row.names = F)















