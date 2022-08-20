## runBestSubsets.R
# finds the subsets of SBC genes that best predict NASH and steatosis

library(leaps)
library(RNOmni)
library(ggplot2)
library(edgeR)
library(bestglm)
library(precrec)

prepDataframe <- function(phenotype, gene_ID, gene_symbol, expr0, covSub) {
  # merge on sample ID
  expr_sub <- expr0[,gene_ID]
  ov = merge(expr_sub, covSub, by = 0, all = F)

  ov = na.omit(ov)
  phe = ov[,paste0(phenotype, "_grp")]

  # get expression data with rows matched to sample ID
  expr = ov[,2:(1+length(gene_ID))]
  colnames(expr) <- gene_symbol

  # INT the expression data
  exprINT <- data.frame(apply(expr, 2, RankNorm))

  return(list(ov = ov, phe = phe, exprINT = exprINT))
}


# function to run best subsets, correct for covariates first
runLeaps <- function(phenotype, gene_ID, gene_symbol, expr0, covSub, precorrect, forcedidx) {
  prep = prepDataframe(phenotype, gene_ID, gene_symbol, expr0, covSub)
  ov = prep$ov
  phe = prep$phe
  exprINT = prep$exprINT

  if(precorrect == T) {
    # pre-correct the counts for covariates
    predictors = apply(exprINT, 2, function(gene) residuals(lm(
          gene ~ ov$Age + ov$isMale + ov$RIN +
          ov$UNIQ_MAP_PERCENT + ov$PCT_INTRONIC_BASES + 
          ov$MEDIAN_3PRIME_BIAS)))
    predictors = as.data.frame(predictors)
  } else {
    # force covariates into the model
    predictors = ov 
    # put genes and covariates together
    predictors[gene_ID] = exprINT
    # remove unneeded columns
    rownames(predictors) = predictors$Row.names 
    predictors[paste0(phenotype, '_grp')] <- NULL
    predictors['Row.names'] <- NULL
    # clean up gene names
    colnames(predictors)[1:length(gene_symbol)] <- gene_symbol
  }

  # run leaps (linear regression)
  bs <- regsubsets(phe ~ ., data = predictors, nbest = 1, 
                    nvmax = ncol(predictors), force.in = forcedidx)

  if (length(forcedidx) == 0) {
    # run leaps (logistic regression)
    # prep input dataframe
    pred_out = data.frame(cbind(predictors, phe))
    colnames(pred_out)[ncol(pred_out)] = pheno
    # run leaps
    bs_logistic = bestglm(pred_out, family = binomial())

    # get AUC from logistic regression
    log_predict = predict(bs_logistic$BestModel, pred_out, type = "response")
    sscurve = evalmod(scores = log_predict, labels = phe)
    aucs = auc(sscurve)
    auc = aucs[aucs$curvetypes == "ROC",]$aucs

    s = data.frame(bs_logistic$Subsets)
    ingenes = s[grepl('*', rownames(s), fixed = T),]
    bic = ingenes[,'BIC']
    ingenes = names(ingenes)[which(ingenes == T)]
    logistic_result = data.frame(phenotype = phenotype,
                                precorrect = precorrect,
                                forcedin = length(forcedidx) > 0,
                                var_in_model = paste(ingenes, collapse = ';'),
                                ngenes_in_model = length(ingenes) - 1 - length(forcedidx),
                                auc = auc,
                                bic = bic)
  }

  return(list(linear = bs, logistic = logistic_result))
}

leapsToDataframe <- function(leapsobj) {
  s <- summary(leapsobj)
  df <- data.frame(s$outmat)
  df["r2"] <- s$rsq
  df["adjr2"] <- s$adjr2
  df["cp"] <- s$cp
  df["bic"] <- s$bic
  df["rss"] <- s$rss
  rownames(df) <- seq(1, nrow(df))
  return(df)
}

getBestModel <- function(leapsDF, phenotype, precorrect, forcedidx) {
  bestrow = leapsDF[which.min(leapsDF$bic),]
  inmodel = names(bestrow)[which(bestrow == "*")]
  n_genes = sum(!(inmodel %in% c("Age", "isMaleTRUE", "RIN", "UNIQ_MAP_PERCENT", 
                                  "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")))
  r2 = bestrow[,"r2"]
  adjr2 = bestrow[,"adjr2"]
  bic = bestrow[,"bic"]

  return(data.frame(phenotype = phenotype,
                    precorrect = precorrect,
                    forcedin = length(forcedidx) > 0,
                    var_in_model = paste(inmodel, collapse = ';'), 
                    ngenes_in_model = n_genes,
                    r2 = r2, 
                    adjr2 = adjr2, 
                    bic = bic))
}

# read in SBC DE results
sbc = read.table("../../data/metTable_subset_serumBiomarkers.txt", sep = "\t", header = T)

# read in adipose expression data
expr0 = read.table("/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt", header = TRUE, row.names = 1)
dge0 = DGEList(expr0)

# restrict to same expressed genes from DE
keepGeneList = read.table("../../data/allExprGenes_adipose.txt")$V1
dge0 = dge0[keepGeneList, , keep.lib.sizes = F]

# calculate CPMs
dge0 = calcNormFactors(dge0)
cpm = data.frame(t(cpm(dge0, log = F)))

# read in phenotype/covariate data
cov = read.csv("../../data/cov.csv", row.names = 1)
cov["IID"] = rownames(cov)

# run best subsets for all histology phenotypes
for (pheno in c("steatosis", "fibrosis", "diagnosis")) {
  leaps_result = data.frame()
  leaps_result_logistic = data.frame()
  # get phenotype status and covariates
  descols = c(paste0(pheno, "_grp"), "IID", "Age", "isMale", "RIN", "UNIQ_MAP_PERCENT", "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")
  covSub = cov[,descols]#, "PRS_BMI")]
  rownames(covSub) = covSub$IID
  covSub = subset(covSub, select = -c(IID))

  # select SBCs DE for this phenotype, the input list of predictor variables
  sbcSub = sbc[sbc$phenotype == pheno,]
  # if steatosis, select only genes that could be used for early detection
  if (pheno == "steatosis") {
    sbcSub = sbcSub[!(sbcSub$gene_ID %in% sbc[sbc$phenotype == "diagnosis",]$gene_ID),]
  }
  
  # run best subsets (pre-correct for covariates)
  leaps_out = runLeaps(pheno, sbcSub$gene_ID, sbcSub$gene_symbol, cpm, covSub, precorrect = T, forcedidx = c())
  leaps_df = leapsToDataframe(leaps_out$linear)
  write.table(leaps_df, file = paste0("../../data/bestSubsets_", pheno, ".txt"), quote = F, row.names = F, sep = '\t')

  leaps_result = rbind(leaps_result, getBestModel(leaps_df, pheno, precorrect = T, forcedidx = c()))
  leaps_result_logistic = rbind(leaps_result_logistic, leaps_out$logistic)

  # run best subsets (put covariates in model and force them in)
  nGenes = length(sbcSub$gene_ID)
  nCov = dim(covSub)[2] - 1
  fi = seq(nGenes + 1, (nGenes + nCov))
  leaps_out_covforced = runLeaps(pheno, sbcSub$gene_ID, sbcSub$gene_symbol, cpm, covSub, precorrect = F, fi)
  leaps_df_covforced = leapsToDataframe(leaps_out_covforced$linear)
  write.table(leaps_df_covforced, file = paste0("../../data/bestSubsets_", pheno, "_covForcedIn.txt"), quote = F, row.names = F, sep = '\t')

  leaps_result = rbind(leaps_result, getBestModel(leaps_df_covforced, pheno, precorrect = F, forcedidx = fi))
  # leaps_result_logistic = rbind(leaps_result_logistic, leaps_out_covforced$logistic)

  write.table(leaps_result, paste0("../../data/", pheno, "_bestmodels.txt"), row.names = F)
  write.table(leaps_result_logistic, paste0("../../data/", pheno, "_bestmodels_logistic.txt"), row.names = F)
}








