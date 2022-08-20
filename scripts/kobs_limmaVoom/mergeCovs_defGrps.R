## mergeCovs_defGrps.R
## import raw data, merge covariates into one table, define groups

# import phenotype data
pheno <- read.csv("../../data/rawpheno/KOBSBaselinePheno.csv", header = TRUE, row.names = 1, na.strings = "#NULL!")
rownames(pheno) <- paste("X", rownames(pheno), sep = "")
# import tech factor data
tf <- read.table("/u/project/pajukant/dzpan29/Juno/KOBS/Complete_RNA_tech_metrics_KOBS_noMT_BaselineAll.txt", header = TRUE, row.names = 1)
rownames(tf) <- paste("X", rownames(tf), sep = "")

# clean up phenotype data
prownames <- rownames(pheno)
pcolnames <- c("sex", "dm_1=yes", "cholesterolmed_preoper", "BMI", 
           "Alaninen_aminotransferase", "Cholesterol", "HDL_cholesterol", "LDL_cholesterol", 
           "Triglycerides", "Fasting_glucose", "Fasting_insulin", "fasting_free_fatty_acids", 
           "Age", "Steatosisgradebaseline", "Fibrosisstagebaseline", "Lobularinflammationbaseline", 
           "Ballooningbaseline", "Diagnosisbaseline", "fenotype_group")
pheno <- data.frame(apply(pheno, 2, function(x) as.numeric(as.character(x))))
rownames(pheno) <- prownames
colnames(pheno) <- pcolnames
# remove "fenotype_group" column
pheno <- pheno[-19]

# define groups
# unhealthy/healthy (defined by sum of all liver histology)
liv_hist <- c("Steatosisgradebaseline", "Fibrosisstagebaseline", 
              #"Lobularinflammationbaseline", "Ballooningbaseline", 
              "Diagnosisbaseline")
pheno$unhealthy_grp <- rowSums(pheno[,liv_hist]) > 0

# function to define groups based on histology
defgrp <- function(hist, unhealthy) {
	grp = hist > 0
	grp[hist == F & unhealthy == T] <- NA
	return(grp)
}

# NASH diagnosis/no NASH diagnosis
pheno$diagnosis_grp <- defgrp(pheno[,"Diagnosisbaseline"], pheno$unhealthy_grp)
# fibrosis/no fibrosis
pheno$fibrosis_grp <- defgrp(pheno[,"Fibrosisstagebaseline"], pheno$unhealthy_grp)
# steatosis/no steatosis
pheno$steatosis_grp <- defgrp(pheno[,"Steatosisgradebaseline"], pheno$unhealthy_grp)
# # ballooning/no ballooning
# pheno$ballooning_grp <- defgrp(pheno[,"Ballooningbaseline"], pheno$unhealthy_grp)
# # inflammation/no inflammation
# pheno$inflammation_grp <- defgrp(pheno[,"Lobularinflammationbaseline"], pheno$unhealthy_grp)

# keep track of how many samples this restricts us to
cc = as.data.frame(apply(pheno[,19:22], 2, function(x) table(x, useNA = "ifany")))
write.table(cc, "../../data/casecontrolratio_histgroups.txt", quote = F, sep = '\t')

# create one merged covariate table (exclude redundant age and sex columns)
pheno_sub <- subset(pheno, select = -c(sex, Age))
cov <- merge(pheno_sub, tf, by = 0)
rownames(cov) <- cov$Row.names
cov <- cov[-1]

# rename sex column
colnames(cov)[36] <- "isMale"
cov$isMale <- cov$isMale == 1

# remove low-quality samples (116, 293, and 328)
cov <- cov[!(rownames(cov) %in% c("X116", "X293", "X328")), ]
pheno0 <- pheno
pheno <- pheno[rownames(pheno) %in% rownames(cov), ]

# # count the number of people with each grade of steatosis
# stect <- data.frame(table(cov$Steatosisgradebaseline))
# colnames(stect) <- c("steatosis_grade", "count")
# write.table(stect, "../../data/steatosis_grade_table_baseline.txt", quote = FALSE, row.names = FALSE)

# make an alternate covariate table for the liver QC metrics (use the same phenotypes per sample)
# import picard metrics
liverpic <- read.table("../../data_liver/picardRNAmetrics_merged_kobs_liver_noMT.txt", header = T, sep = '\t')
liverpic["sample"] <- paste("X", liverpic$SAMPLE, sep = "")
# import RIN data
liverRIN <- read.csv("/u/project/pajukant/malvarez/KOBS_liver_align/data/raw/KOBS_adipose_and_liver_data/KOBS_liver_plate.csv")
liverRIN <- liverRIN[,c("sample.ID", "RINe")]
colnames(liverRIN) <- c("sample", "RIN")
liverRIN["sample"] <- paste("X", as.character(liverRIN$sample), sep = "")
# import STAR QC metrics
liverSTAR <- read.table("/u/project/pajukant/malvarez/KOBS_liver_align/data/processed/qc/map_stats/map_stats.subject.txt", header = T)
liverSTAR <- liverSTAR[,c("sample", "pct_uniq_map")]
colnames(liverSTAR)[2] <- "UNIQ_MAP_PERCENT"
liverSTAR["sample"] <- paste("X", liverSTAR$sample, sep = "")
# crop the adipose covariates to just the ones in common with liver
liverCov <- pheno0
# rename sex column
colnames(liverCov)[1] <- "isMale"
liverCov$isMale <- liverCov$isMale == 1
liverCov <- liverCov[,c("Age", "isMale", "steatosis_grp", "fibrosis_grp", 
              # "inflammation_grp", "ballooning_grp", 
              "diagnosis_grp",
              "unhealthy_grp")]
liverCov["sample"] <- rownames(liverCov)
# paste RIN onto this table
liverCov <- merge(liverCov, liverRIN, by = "sample")
# paste picard tech factors onto this table
liverCov <- merge(liverCov, liverpic, by = "sample")
# paste STAR QC metrics onto this table
liverCov <- merge(liverCov, liverSTAR, by = "sample", how = "left")

# write files out
write.csv(pheno, file = "../../data/pheno.csv")
write.csv(cov, file = "../../data/cov.csv")
write.csv(liverCov, file = "../../data_liver/cov_liver.csv")






