## associateCCDC80.R
# associate CCDC80 adipose expression in KOBS with traits we think might be related

library(ggplot2)
library(RNOmni)

setwd("~/src/adiDys/bin/scripts/")

# import count and covariate data
countfile = "../../data/rawcounts/BaselineAll_KOBS_geneCount_noMT.txt"
counts <- read.table(countfile, header = TRUE, row.names = 1)
nsamples <- ncol(counts)
cov <- read.csv("../../data/cov.csv", header = TRUE, row.names = 1)

# select traits we're interested in
cov <- cov[,c("dm_1.yes", "Triglycerides", "Age", "isMale", "RIN", 
              "UNIQ_MAP_PERCENT", "PCT_INTRONIC_BASES", "MEDIAN_3PRIME_BIAS")]

# select SBC genes
sbc <- read.table("../../data/metTable_subset_serumBiomarkers.txt", sep = '\t', header = T)[,c("gene_ID", "gene_symbol")]
sbc <- unique(sbc)
counts <- counts[rownames(counts) %in% sbc$gene_ID,]
counts["gene_ID"] <- rownames(counts)
counts <- merge(counts, sbc, by = "gene_ID")
rownames(counts) <- counts$gene_symbol
counts <- counts[2:(2+nsamples-1)]

# merge the dataframes to put the samples in the same order
cov["sample"] <- rownames(cov)
counts <- data.frame(t(counts))
counts["sample"] <- rownames(counts)
m <- merge(cov, counts, by = "sample", type = "left")
m["dm_1.yes"] = as.factor(m$dm_1.yes)
m <- na.omit(m)

# INT expression data and correct for covariates
expr <- m[,colnames(m) %in% sbc$gene_symbol]
expr_INT <- data.frame(apply(expr, 2, RankNorm))
expr_INT_cov <- data.frame(apply(expr_INT, 2, function(gene) residuals(lm(
  gene ~ m$Age + m$isMale + m$RIN + m$UNIQ_MAP_PERCENT + 
    m$PCT_INTRONIC_BASES + m$MEDIAN_3PRIME_BIAS))))
m <- cbind(m[1:9], expr_INT_cov)

# INT TG data and correct for age and sex
tg_INT <- RankNorm(m$Triglycerides)
tg_INT_cov <- residuals(lm(tg_INT ~ m$Age + m$isMale))
m["TG_INT_COV"] <- tg_INT_cov

# associate SBC genes with corrected TG values
tg_assoc <- data.frame(matrix(nrow = 3, ncol = 0))
for (gene in sbc$gene_symbol) {
  genevec <- m[[gene]]
  fit <- lm(data = m, TG_INT_COV ~ genevec)
  coef <- coefficients(summary(fit))
  pval <- coef["genevec", "Pr(>|t|)"]
  slope <- coef["genevec", "Estimate"]
  tg_assoc[gene] <- c(gene, pval, slope)
}
rownames(tg_assoc) <- c("gene", "p", "slope")
tg_assoc <- data.frame(t(tg_assoc))
tg_assoc["phenotype"] <- "TG"

# associate SBC genes with T2D diagnoses
t2d_assoc <- data.frame(matrix(nrow = 3, ncol = 0))
for (gene in sbc$gene_symbol) {
  genevec <- m[[gene]]
  fit <- glm(data = m, dm_1.yes ~ genevec, family = "binomial")
  coef <- coefficients(summary(fit))
  pval <- coef["genevec", "Pr(>|z|)"]
  slope <- coef["genevec", "Estimate"]
  t2d_assoc[gene] <- c(gene, pval, slope)
}
rownames(t2d_assoc) <- c("gene", "p", "slope")
t2d_assoc <- data.frame(t(t2d_assoc))
t2d_assoc["phenotype"] <- "T2D"

write.table(tg_assoc, file = "../../data/SBC_TG_association.txt", sep = '\t', quote = F, row.names = F)
write.table(t2d_assoc, file = "../../data/SBC_T2D_association.txt", sep = '\t', quote = F, row.names = F)


# # associate CCDC80 expression with LPL expression
# 
# # CCDC80 = ENSG00000091986
# # LPL = ENSG00000175445
# ggplot(data = m, aes(x = ENSG00000091986, y = ENSG00000175445)) + 
#   geom_point() +
#   xlab("CCDC80") +
#   ylab("LPL")
# 
# fitLPL <- lm(data = m, ENSG00000175445 ~ ENSG00000091986)
# summary(fitLPL)
