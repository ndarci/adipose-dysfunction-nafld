# run a permutation test testing for significance of best subsets models
library(RNOmni)
library(edgeR)

# define number of permutations we'll do
B = 100000

cleanupDataFrame <- function(phenotype, expr0, covSub, permgenelist, snum) {
    # get some random genes for the permutation
    randgenes = sample(permgenelist, snum)
    # set up these genes for a linear regression
    if (all(randgenes %in% names(expr0))) {
    	expr_sub <- data.frame(expr0[,randgenes])
    } else {
    	return('failed')
    }
    rownames(expr_sub) = rownames(expr0)
    ov = merge(expr_sub, covSub, by = 0, all = F)
    # remove NA and select phenotype data
    ov = na.omit(ov)
    phe = ov[,paste0(phenotype, "_grp")]
    return(list(ov = ov, phe = phe, randgenes = randgenes))
}

# general fx to run permutation test
permutationTest <- function(phenotype, bestmodel, expr0, covSub, permgenelist) {
  # remember number of genes in best model
  snum = bestmodel[,'ngenes_in_model']
  # remember r2 of each permuted model
  r2vec = rep(NA, B)
  # iterate over all permutations
  for (i in seq(1, B)) {
  	while(TRUE) {
	  	clean = cleanupDataFrame(phenotype, expr0, covSub, permgenelist, snum)
	  	if(length(clean) > 1) break
  	}
  	ov = clean$ov 
  	phe = clean$phe
  	randgenes = clean$randgenes
    # get expression data with rows matched to sample ID
    expr = data.frame(ov[,2:(1+length(randgenes))])
    colnames(expr) <- randgenes
    # INT the expression data
    exprINT <- data.frame(apply(expr, 2, RankNorm))
    # regress out covariates
    exprINT_cov <- apply(exprINT, 2, function(gene) residuals(lm(
            gene ~ ov$Age + ov$isMale + ov$RIN +
            ov$UNIQ_MAP_PERCENT + ov$PCT_INTRONIC_BASES + 
            ov$MEDIAN_3PRIME_BIAS)))
    exprINT_cov = as.data.frame(exprINT_cov)
    fit = lm(phe ~ ., data = exprINT_cov)
    r2lm = summary(fit)$r.squared
    r2vec[i] = r2lm
  }
  # get proportion of permuted r2 values more extreme than real r2
  realr2 = bestmodel[,"r2"]
  pval = sum(r2vec > realr2) / length(r2vec)
  return(pval)
}

permutationTest_cov <- function(phenotype, bestmodel, expr0, covSub, permgenelist, covOut) {
  # remember number of genes in best model
  snum = bestmodel[,'ngenes_in_model']
  # remember r2 of each permuted model
	r2vec = rep(NA, B)
	# iterate over all permutations
	for (i in seq(1, B)) {
  	while(TRUE) {
	  	clean = cleanupDataFrame(phenotype, expr0, covSub, permgenelist, snum)
	  	if(length(clean) > 1) break
  	}
		ov = clean$ov 
		phe = clean$phe
  	randgenes = clean$randgenes
		# get expression and phenotype data with rows matched to sample ID
		expr = data.frame(ov[,2:(1+length(randgenes))])
		colnames(expr) <- randgenes

		# INT the expression data
		exprINT <- data.frame(apply(expr, 2, RankNorm))
		# update expression data in ov table
		ov[randgenes] = exprINT 

		# add covariates to the model
		rownames(ov) = ov$Row.names
		ov[paste0(phenotype, "_grp")] <- NULL
		ov["Row.names"] <- NULL

		# remove covariates leaps didn't choose
		ov = ov[,!(names(ov) %in% covOut)]

		# fit linear model
		fit = lm(phe ~ ., data = ov)
		r2lm = summary(fit)$r.squared
		r2vec[i] = r2lm
	}
	# get proportion of permuted r2 values more extreme than real r2
	realr2 = bestmodel[,"r2"]
	pval = sum(r2vec > realr2) / length(r2vec)
	return(pval)
}

# getCovOutList <- function(bestmodel) {
# 	result = names(bestmodel)[which(bestmodel == " ")]
# 	return(result)
# }

# read in adipose expression data
expr0 = read.table("/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt", header = TRUE, row.names = 1)
dge0 = DGEList(expr0)

# restrict to same expressed genes from DE
permgenelist = read.table("../../data/allExprGenes_adipose.txt")$V1
dge0 = dge0[permgenelist, , keep.lib.sizes = F]

# calculate CPMs
dge0 = calcNormFactors(dge0)
expr0 = data.frame(t(cpm(dge0, log = F)))

# read in covariate data
cov = read.csv("../../data/cov.csv", row.names = 1)
cov["IID"] = rownames(cov)


result = data.frame()
for (pheno in c("steatosis", "fibrosis", "diagnosis")) {
	# select covariates and target phenotype
	covSub = cov[,c("IID", paste0(pheno, "_grp"), "Age", "isMale", 
	                "RIN", "UNIQ_MAP_PERCENT", "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")]
	covSub = subset(covSub, select = -c(IID))

	# read in best subsets results
	leaps = read.table(paste0("../../data/", pheno, "_bestmodels.txt"), header = T)

	# run permutation test for pre-corrected model
	pval = permutationTest(pheno, leaps[grepl("pre-corrected", leaps$name),], expr0, covSub, permgenelist)

	# run permutation test for forced in model
	pval_forced = permutationTest_cov(pheno, leaps[grepl("forced in", leaps$name),], expr0, covSub, permgenelist, c())

	result = rbind(result, data.frame(phenotype = pheno,
																		pval_precorrected = pval,
																		pval_covforced = pval_forced))
}
write.table(result, "../../data/bestSubsets_permutationTestResult.txt", quote = F, row.names = F)





