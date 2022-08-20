## runLimmaVoom.R
## run limma-voom for DE analysis

tissue <- commandArgs(trailingOnly = T)[1]

library(edgeR)
library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)

#setwd('~/src/adiDys/bin/scripts')

# import count and covariate data (mitochondrial reads already removed)
if (tissue == "adipose") {
  countfile = "/u/project/pajukant/dzpan29/Juno/KOBS/BaselineAll_KOBS_geneCount_noMT.txt"
  counts <- read.table(countfile, header = TRUE, row.names = 1)
  cov <- read.csv("../../data/cov.csv", header = TRUE, row.names = 1)
  
  cov["dm_1.yes"] = as.logical(as.numeric(cov$dm_1.yes))
} else { # tissue is liver
  countfile = "../../data_liver/gene_counts.match.txt"
  counts <- read.table(countfile, header = TRUE, row.names = 1)
  counts <- counts[,7:ncol(counts)]
  cov <- read.csv("../../data_liver/cov_liver.csv", header = TRUE, row.names = 2)
}

# # demonstrate how much highly expressed genes can drive results
# totReads <- sum(rowSums(counts))
# top500readcount <- sum(tail(sort(rowSums(counts)), n = 500))
# # what proportion of the total reads are mapped to the top 500 most expressed genes? (~58k total genes)
# top500readcount / totReads

# # check the percent of male individuals in the cohort
# tsex <- table(cov$isMale)
# tsex
# tsex[[2]] / (tsex[[1]] + tsex[[2]])
# # check the sample sizes by group and sex for steatosis
# table(cov[,c("steatosis_grp", "isMale")])

# create DGE object and filter out lowly expressed genes
dge0 <- DGEList(counts)
# keep only those genes with nonzero expression in at least 90% of the samples
keepGenes <- rowSums(dge0$counts > 0) >= 0.9 * length(colnames(dge0$counts))
# how many genes pass the filter?
table(keepGenes)
# update the DGE object with only genes we want
dge0 <- dge0[keepGenes, , keep.lib.sizes = FALSE]

stripVersion <- function(idlist) {
        strid = strsplit(idlist, ".", fixed = T)
        stripped = sapply(strid, "[[", 1)
        return(stripped)
}

# save the list of all expressed genes
exprgenes = rownames(dge0$counts)
if(tissue == "liver") {
	exprgenes = stripVersion(exprgenes)
}
write.table(exprgenes, file = paste0("../../data/allExprGenes_", tissue, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# run several independent DE experiments with different groupings of the samples
# sex-specific test needs a different design matrix WITHOUT sex included as a covariate
getDesign <- function(covar, grp, label, type) {
  # determine appropriate covariates and group comparisons based on design type
  if(type == "hist") {
      design <- model.matrix(~ grp + Age + isMale + RIN + UNIQ_MAP_PERCENT + PCT_INTRONIC_BASES + MEDIAN_3PRIME_BIAS, data = covar)
      colnames(design)[2] <- paste("yes", label, sep = "")
      } else if(type == "sex") {
      design <- model.matrix(~ grp + Age + RIN + UNIQ_MAP_PERCENT + PCT_INTRONIC_BASES + MEDIAN_3PRIME_BIAS, data = covar)
      colnames(design)[2] <- paste("yes", label, sep = "")
      } else { return("failed") }
  return(design)
}

options(na.action = 'na.pass')
runDE <- function(dge, design) {
  # filter out missing data
  design <- na.omit(design)
  keepSamples <- rownames(design)
  keepSamples <- intersect(keepSamples, rownames(dge$samples))
  dge <- dge[,keepSamples]
  design <- design[keepSamples,]
  # normalize by library size
  dge <- calcNormFactors(dge)
  
  # # save the TMM-normalized counts
  # test <- colnames(design)[2]
  # test <- substr(test, 4, nchar(test))
  # if (tissue == "adipose") {
  #   tmmfilename <- paste("../../data/TMMnormCPM_", test, "_", tissue, ".txt", sep = "")
  # } else {
  #   tmmfilename <- paste("../../data_liver/TMMnormCPM_", test, "_", tissue, sep = "")
  # }
  # write.table(cpm(dge, normalized.lib.sizes = T), file = tmmfilename, quote = F)
  
  # use voom to log2 transform the counts
  v <- voom(dge, design)
  # create a linear model representing the effects of each covariate and the experimental group
  fit <- lmFit(v, design)
  # compute Bayes statistics
  fit2 <- eBayes(fit)
  # get top DE genes
  # topGenes <- topTable(fit2, coef = 2, n = Inf, p.value = 0.05, sort.by = "P")
  topGenes <- topTable(fit2, coef = 2, n = Inf, p.value = Inf, sort.by = "P")
  return(topGenes)
  # return(list(topGenes = topGenes, fit = fit2))
}

# generate design matrices for histology traits
fibDes <- getDesign(cov, cov$fibrosis_grp, "Fibrosis", "hist")
steDes <- getDesign(cov, cov$steatosis_grp, "Steatosis", "hist")
diaDes <- getDesign(cov, cov$diagnosis_grp, "NASHdiagnosis", "hist")
# balDes <- getDesign(cov, cov$ballooning_grp, "Ballooning", "hist")
# infDes <- getDesign(cov, cov$inflammation_grp, "Inflammation", "hist")
# unhDes <- getDesign(cov, cov$unhealthy_grp, "Unhealthy", "hist")

#t2dDes <- getDesign(cov, cov$dm_1.yes, "T2D", "hist")

# run DE analysis for histology traits
fibDE <- runDE(dge0, fibDes)
steDE <- runDE(dge0, steDes)
diaDE <- runDE(dge0, diaDes)
# balDE <- runDE(dge0, balDes)
# infDE <- runDE(dge0, infDes)
# unhDE <- runDE(dge0, unhDes)

#t2dDE <- runDE(dge0, t2dDes)

# # run DE analysis separately for males and females (just steatosis for now)
# cov_male <- cov[cov$isMale, ]
# cov_female <- cov[!cov$isMale, ]
# steDes_male <- getDesign(cov_male, cov_male$steatosis_grp, "Steatosis", "sex")
# steDes_female <- getDesign(cov_female, cov_female$steatosis_grp, "Steatosis", "sex")
# steDE_male <- runDE(dge0, steDes_male)
# steDE_female <- runDE(dge0, steDes_female)
# 
# # run DE analysis comparing males and females
# sexDes <- getDesign(cov, cov$isMale, "Male", "sex")
# sexDE <- runDE(dge0, sexDes)

# # run sex-specific DE again, using only genes that were DE in STEATOSIS and SEX (first pass)
# # this time use the INTERACTION term
# steSexGenes <- rownames(sexDE)[rownames(sexDE) %in% rownames(steDE)]
# dge_ss <- dge0[steSexGenes, , keep.lib.sizes = FALSE]
# dge_ss <- calcNormFactors(dge_ss)
# 
# # glm(steatosis ~ sex*G , data = mydataframe , family = "binomial")
# 
# # steSexDes <- getDesign(cov, cov$isMale, "Male", "sex_interaction")
# # steSexDE <- runDE(dge_ss, steSexDes)

# write DE gene lists and topTables to files
writeDE <- function(tt, label) {
  if (tissue == "adipose") {
    DElistFilename = paste("../../data/DEgenes_", label, "_", tissue, ".txt", sep = "")
    TTfilename = paste("../../data/topTable_", label, "_", tissue, ".txt", sep = "")
  } else {
    DElistFilename = paste("../../data_liver/DEgenes_", label, "_", tissue, ".txt", sep = "")
    TTfilename = paste("../../data_liver/topTable_", label, "_", tissue, ".txt", sep = "")
  }
  # write the stats for all results
  nomfn = paste0(substr(TTfilename, 1, nchar(TTfilename)-4), '_nominal.txt')
  write.table(tt, file = nomfn, row.names = TRUE, col.names = TRUE, quote = TRUE)
  # write the stats for significant results 
  tt = tt[tt$adj.P.Val < 0.05,]
  write.table(rownames(tt), file = DElistFilename, row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(tt, file = TTfilename, row.names = TRUE, col.names = TRUE, quote = TRUE)
}
writeDE(fibDE, "fibrosis")
writeDE(steDE, "steatosis")
writeDE(diaDE, "diagnosis")
# writeDE(balDE, "ballooning")
# writeDE(infDE, "inflammation")
# writeDE(unhDE, "unhealthy")

if (tissue == "liver") {
  # do a sex-stratified DE for steatosis in liver
  cov_m = cov[cov$isMale == T,]
  cov_f = cov[cov$isMale == F,]
  print(dim(cov_m))
  print(dim(cov_f))
  steDes_m = getDesign(cov_m, cov_m$steatosis_grp, "Steatosis", "sex")
  steDes_f = getDesign(cov_f, cov_f$steatosis_grp, "Steatosis", "sex")
  steDE_m = runDE(dge0, steDes_m)
  steDE_f = runDE(dge0, steDes_f)
  writeDE(steDE_m, "steatosis_male")
  writeDE(steDE_f, "steatosis_female")
}

#writeDE(t2dDE, "T2D")

# writeDE(steDE_male, "steatosis_male")
# writeDE(steDE_female, "steatosis_female")
# writeDE(sexDE, "sex")
# writeDE(steSexDE, "steatosis_sex_interaction")

# # check overlap between DE genes in steatosis_female and sex comparison
# s <- rownames(sexDE)
# sf <- rownames(steDE_female)
# ov <- s[s %in% sf]

# # double-check direction of logFC by plotting expression
# coolGene <- "ENSG00000170525"
# steatSamps <- rownames(na.omit(cov[cov$steatosis_grp == TRUE, ]))
# noSteatSamps <- rownames(na.omit(cov[cov$steatosis_grp == FALSE, ]))
# infSamps <- rownames(na.omit(cov[cov$inflammation_grp == TRUE, ]))
# noInfSamps <- rownames(na.omit(cov[cov$inflammation_grp == FALSE, ]))
# steatCPM <- cpms[coolGene, steatSamps]
# noSteatCPM <- cpms[coolGene, noSteatSamps]
# infCPM <- cpms[coolGene, infSamps]
# noInfCPM <- cpms[coolGene, noInfSamps]
# df <- data.frame(grp = c(rep("Steatosis", length(steatCPM)), rep("No steatosis", length(noSteatCPM)),
#                          rep("Inflammation", length(infCPM)), rep("No inflammation", length(noInfCPM))),
#                  expr = c(steatCPM, noSteatCPM, infCPM, noInfCPM))
# lev_order <- factor(df$grp, levels = c("Steatosis", "No steatosis", "Inflammation", "No inflammation"))
# ggplot(data = df, aes(x = lev_order, y = expr)) + geom_boxplot() + scale_y_continuous(trans = "log2") + ylab("CPM (log2 scale)") + xlab("group") + ggtitle("Expression of PFKFB3")






