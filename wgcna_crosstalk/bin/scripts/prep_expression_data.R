tissue = commandArgs(trailingOnly = T)[1]

library(edgeR)
library(RNOmni)

# set global data.frame options
options(stringsAsFactors = FALSE)

# read in expression data and covariates/phenotypes
if (tissue == "adipose") {
  	countfile = "/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt"
  	counts <- read.table(countfile, header = TRUE, row.names = 1)
	cov <- read.csv("/u/project/pajukant/nikodm/kobs_limmaVoom/data/cov.csv", header = TRUE, row.names = 1)
	cov["dm_1.yes"] = as.logical(as.numeric(cov$dm_1.yes))
} else { # tissue is liver
	countfile = "/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/gene_counts.match.txt"
	counts <- read.table(countfile, header = TRUE, row.names = 1)
	counts <- counts[,7:ncol(counts)]
	cov <- read.csv("/u/project/pajukant/nikodm/kobs_limmaVoom/data_liver/cov_liver.csv", header = TRUE, row.names = 2)
}

# create DGE object
y0 <- DGEList(counts)

# keep only those genes with nonzero counts in at least 90% of the samples
keepGenes <- rowSums(y0$counts > 0) >= 0.9 * length(colnames(y0$counts))
# how many genes pass the filter?
table(keepGenes)
y0 <- y0[keepGenes, , keep.lib.sizes = FALSE]

# define covariates we care about
covlist = c("Age", "isMale", "RIN", "UNIQ_MAP_PERCENT", "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")

# restrict data to invidivials with data for each of these covariates, and across liver/adipose
validsamples = read.table("../../data/KOBS_adipose_liver_overlap_samples_RNAseq.txt")$V1
cov_clean = na.omit(cov[validsamples, covlist])
y0 <- y0[, validsamples, keep.lib.sizes = FALSE]

# get norm factors to calculate CPM
y0 <- calcNormFactors(y0)

# calculate CPMs
cpm = cpm(y0, log = F)

# INT the CPMs to meet linear regression normality assumption
cpm_INT = data.frame(apply(cpm, 1, RankNorm))

# regress out tech factors
covlist_formula = paste0("cov_clean$", covlist)
cpm_INT_resid = data.frame(apply(cpm_INT, 2, function(gene) 
				residuals(
				lm(
				as.formula(
				paste0("gene~", paste(covlist_formula, collapse = "+"))
				)))))

# # do a log transform
# addme = abs(min(cpm_resid))+1
# cpm_resid_log = data.frame(apply(cpm_resid, 2, function(gene) log2(gene+addme)))

# # plot distribution of some random genes to check normality
# library(ggplot2)
# set.seed(22222)
# randgenes = sample(names(cpm_resid), 10)
# for (gene in randgenes) {
# 	p = ggplot(mapping = aes(x = cpm_resid_log[,gene])) + geom_histogram(bins = 50) + xlab(gene) + ylab("Count")
# 	ggsave(paste0('../../fig/norm_check_log_', gene, '.png'), p)
# }

# INT the corrected, filtered, CPMs to meet Pearson correlation normality assumption
# cpm_INT_resid_INT = cpm_INT_resid
cpm_INT_resid_INT = data.frame(apply(cpm_INT_resid, 2, RankNorm))

# use the covariates that are constant per individual and don't depend on tissue (no tech factors)
cov_individual = read.csv("/u/project/pajukant/nikodm/kobs_limmaVoom/data/cov.csv", header = TRUE, row.names = 1)
cov_individual["dm_1.yes"] = as.logical(as.numeric(cov_individual$dm_1.yes))
phenolist = c("BMI", "dm_1.yes", 
	# "Cholesterol", "HDL_cholesterol", "LDL_cholesterol", 
      "Triglycerides",
      "Fasting_glucose", #"Fasting_insulin", "fasting_free_fatty_acids",
      "Steatosisgradebaseline", "Fibrosisstagebaseline", "Diagnosisbaseline",
      "steatosis_grp", "fibrosis_grp", "diagnosis_grp",
      "cholesterolmed_preoper"
      # "Lobularinflammationbaseline", "Ballooningbaseline",  "unhealthy_grp"
      )
datTraits = merge(cov_individual[validsamples,phenolist], cov_clean, by = 0)
# remove rows with NA in the columns that don't have a ton of missing values
datTraits = datTraits[complete.cases(datTraits[,phenolist[!(phenolist %in% c("steatosis_grp", "fibrosis_grp", "diagnosis_grp"))]]),]
# datTraits = na.omit(datTraits)
rownames(datTraits) = datTraits$Row.names 
datTraits = datTraits[,!(colnames(datTraits) %in% c("Row.names"))]

# correct glucose for T2D
datTraits['Fasting_glucose'] = residuals(lm(Fasting_glucose ~ dm_1.yes, data = datTraits))

# INT the (numerical) phenotypes
phenolist_num = c("Age", "BMI", "Triglycerides", "Fasting_glucose")
datTraits_INT = data.frame(apply(datTraits[,phenolist_num], 2, RankNorm))
datTraits = cbind(datTraits[,!(colnames(datTraits) %in% phenolist_num)], datTraits_INT)

# # add best NAFLD and BMI PRS to the trait data
# prs.bmi = read.table('/u/project/pajukant/nikodm/kobs_prs/data/kobs.0.4.profile', header = T)
# prs.nafld = read.table('/u/project/pajukant/nikodm/kobs_prs/data/kobs.NAFLD_imp.0.4.profile', header = T)
# colnames(prs.bmi)[6] = 'PRS_BMI'
# colnames(prs.nafld)[6] = 'PRS_NAFLD'
# prs.bmi['IID'] = paste0('X', prs.bmi$IID)
# prs.nafld['IID'] = paste0('X', prs.nafld$IID)
# datTraits = merge(datTraits, prs.bmi[,c('IID', 'PRS_BMI')], by.x = 0, by.y = 'IID')
# datTraits = merge(datTraits, prs.nafld[,c('IID', 'PRS_NAFLD')], by.x = 'Row.names', by.y = 'IID')
# rownames(datTraits) = datTraits$Row.names 
# datTraits = datTraits[,!(colnames(datTraits) %in% c("Row.names"))]

# write this clean expression dataframe out along with the covariates
datExpr = cpm_INT_resid_INT[rownames(datTraits),]

datTraits = datTraits[,c(covlist, names(datTraits)[!(names(datTraits) %in% covlist)])]

save(datExpr, datTraits, file = paste0("../../data/", tissue, "_clean_counts_cov.RData"))




